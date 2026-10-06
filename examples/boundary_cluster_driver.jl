#!/usr/bin/env julia
# Cluster → alpha-shape boundary → boundary statistics CLI for SMLMClustering.
#
# Loads a real SMLD `.h5` dataset (via SMLMAnalysis.load_smld), runs DBSCAN
# clustering, extracts alpha-shape cluster boundaries, computes per-cluster
# boundary statistics, prints a summary, and optionally renders cluster
# visualizations — the same pipeline as cluster_driver.jl / boundary_driver.jl,
# but with `smld_path`/`results_dir` as CLI arguments and every other knob
# read from a params.toml file instead of hardcoded literals.
#
# Usage:
#   julia --project=. boundary_cluster_driver.jl \
#       --smld <path.h5> --out <results_dir> [--params <bc_params.toml>]
#
# bc_params.toml (see dev/boundary/bc_params.toml for a full worked example):
#   per_dataset   = false   # shared by [dbscan] and [boundary] unless overridden
#   plotting      = true
#   pixel_size_nm = 1.0
#   [roi]                    # optional — omit to skip ROI filtering
#   x_min = 12.5
#   x_max = 17.5
#   y_min = 12.5
#   y_max = 17.5
#   [dbscan]                 # -> DBSCANConfig(; kw...)
#   eps_nm = 25.0
#   min_points = 3
#   [boundary]                # -> BoundaryClustersConfig(; kw...)
#   shrink = 0.05
#   [plot]                    # -> ClusterPlotConfig(; kw...)
#   colormap = "inferno"
#   color_by = "id"

using TOML
using Dates                  # now
using SMLMAnalysis           # load_smld
using SMLMData                # filter_roi
using SMLMClustering          # cluster, boundary_clusters, boundary_statistics, …
using CairoMakie, SMLMRender  # GaussianRender (for the optional plotting block)
using Statistics               # mean, std

const _USAGE = """
usage: boundary_cluster_driver.jl --smld PATH --out DIR [--params FILE]
"""

# Built-in defaults, used when --params is omitted — mirrors
# cluster_driver.jl's current hardcoded values (no [roi] = no ROI filtering).
const _DEFAULT_PARAMS = Dict{String,Any}(
    "per_dataset"   => false,
    "plotting"      => true,
    "pixel_size_nm" => 1.0,
    "dbscan"        => Dict{String,Any}("eps_nm" => 25.0, "min_points" => 3),
    "boundary"      => Dict{String,Any}("shrink" => 0.05),
    "plot"          => Dict{String,Any}("colormap" => "inferno", "color_by" => "id",
                                         "categorical" => true, "show_centroids" => false),
)

function _parse_args(args)
    out = Dict{String,Any}("smld" => nothing, "out" => nothing, "params" => nothing)
    i = 1
    while i <= length(args)
        a = args[i]
        if a in ("--smld", "--out", "--params")
            i + 1 <= length(args) || error("$a requires a value\n$_USAGE")
            out[a[3:end]] = args[i + 1]; i += 2
        elseif a in ("-h", "--help")
            print(stdout, _USAGE); exit(0)
        else
            error("unknown arg: $a\n$_USAGE")
        end
    end
    for k in ("smld", "out")
        out[k] === nothing && error("--$k is required\n$_USAGE")
    end
    return out
end

# Section => Symbol-keyed kwargs, for splatting into a Base.@kwdef config
# constructor (unknown keys raise a clear error from the constructor itself).
_symbol_kwargs(d::AbstractDict) = Dict{Symbol,Any}(Symbol(k) => v for (k, v) in d)

# Print `name`: mean +/- SEM over the non-NaN entries of `values` (same
# formula as cluster_driver.jl / boundary_driver.jl: std / sqrt(n), i.e. SEM).
function print_stat(name::AbstractString, values::AbstractVector{<:Real},
                     units::AbstractString = "")
    valid = filter(!isnan, values)
    if isempty(valid)
        println(name, ": no valid values")
    else
        n = length(valid)
        println(name, ": ", round(mean(valid), digits = 4),
                " +/- ", round(std(valid) / sqrt(n), digits = 4),
                units == "" ? "" : " " * units, "  (n = ", n, ")")
    end
end

function main()
    opts = _parse_args(ARGS)
    raw  = opts["params"] === nothing ? _DEFAULT_PARAMS : TOML.parsefile(opts["params"])

    per_dataset   = get(raw, "per_dataset", false)
    plotting      = get(raw, "plotting", true)
    pixel_size_nm = get(raw, "pixel_size_nm", 1.0)

    date_time = now()

    @info "loading SMLD" smld_path = opts["smld"]
    smld = load_smld(opts["smld"])

    if haskey(raw, "roi")
        roi = raw["roi"]
        smld = filter_roi(smld, (roi["x_min"], roi["x_max"]), (roi["y_min"], roi["y_max"]))
    end

    # ── Clustering (DBSCAN) ───────────────────────────────────────────────
    dbscan_kw = _symbol_kwargs(get(raw, "dbscan", Dict{String,Any}()))
    haskey(dbscan_kw, :per_dataset) || (dbscan_kw[:per_dataset] = per_dataset)
    cfg_dbscan = DBSCANConfig(; dbscan_kw...)
    (smld_out, info_dbscan) = cluster(smld, cfg_dbscan)
    println(smld_out)
    println(info_dbscan)

    # ── Cluster boundaries (alpha-shape) ──────────────────────────────────
    boundary_kw = _symbol_kwargs(get(raw, "boundary", Dict{String,Any}()))
    haskey(boundary_kw, :per_dataset) || (boundary_kw[:per_dataset] = per_dataset)
    cfg_bnd = BoundaryClustersConfig(; boundary_kw...)
    (_, binfo) = boundary_clusters(smld_out, cfg_bnd)
    println(binfo)
    # binfo.boundaries_x[j] / binfo.boundaries_y[j] trace the boundary polygon
    # of cluster binfo.cluster_keys[j] in μm (closed: first point == last point).

    # ── Cluster boundary statistics ───────────────────────────────────────
    (_, sinfo) = boundary_statistics(smld_out, binfo; per_dataset = per_dataset,
                                      algorithm = info_dbscan.algorithm)
    println()
    println(sinfo)
    println(fieldnames(typeof(sinfo)))

    println()
    println("date_time: ", date_time)
    println("Number of input points: ", sinfo.n_locs_in)
    println("Number of clustered points: ", sinfo.n_clustered)
    println("Number of unclustered points: ", sinfo.n_noise)
    println("Clustering algorithm: ", sinfo.algorithm)
    println("Shrink factor: ", sinfo.shrink)
    println("Number of clusters: ", sinfo.n_clusters)

    print_stat("Area", sinfo.area, "μm²")
    print_stat("Area (convex)", sinfo.areaConvex, "μm²")
    print_stat("Area (tight)", sinfo.areaTight, "μm²")
    print_stat("Perimeter", sinfo.perimeter, "μm")
    print_stat("Perimeter (convex)", sinfo.perimeterConvex, "μm")
    print_stat("Perimeter (tight)", sinfo.perimeterTight, "μm")
    print_stat("Equivalent radius", sinfo.equiv_radius, "μm")
    print_stat("Compactness", sinfo.compactness)
    print_stat("Circularity", sinfo.circularity)
    print_stat("Convexity", sinfo.convexity)
    print_stat("Solidity", sinfo.solidity)
    print_stat("Points per cluster", sinfo.n_pts_per_cluster)
    print_stat("Points per area", sinfo.n_pts_per_area, "1/μm²")
    print_stat("Sigma (intracluster distance std)", sinfo.sigma_actual, "μm")
    print_stat("Cluster width", sinfo.cluster_width, "μm")
    print_stat("Min center-to-center distance", sinfo.min_c2c_dists, "μm")
    println("Min center-to-center distance (overall): ",
            isnan(sinfo.min_c2c_dist) ? "NaN" : round(sinfo.min_c2c_dist, digits = 4), " μm")
    print_stat("Min edge-to-edge distance", sinfo.min_e2e_dists, "μm")
    println("Min edge-to-edge distance (overall): ",
            isnan(sinfo.min_e2e_dist) ? "NaN" : round(sinfo.min_e2e_dist, digits = 4), " μm")
    print_stat("Nearest-neighbor distance (within cluster)",
               reduce(vcat, sinfo.nn_within_clusters; init = Float64[]), "μm")
    println("Elapsed time: ", round(sinfo.elapsed_s * 1e3, digits = 1), " ms")

    if plotting
        results_dir = abspath(opts["out"])
        mkpath(results_dir)
        stem = splitext(basename(opts["smld"]))[1]

        plot_kw = _symbol_kwargs(get(raw, "plot", Dict{String,Any}()))
        for sym in (:colormap, :color_by, :boundary_color)
            haskey(plot_kw, sym) && (plot_kw[sym] = Symbol(plot_kw[sym]))
        end
        plot_kw[:pixel_size_nm] = pixel_size_nm

        # ── Cluster visualization ───────────────────────────────────────────
        png_path = joinpath(results_dir, stem * "_Gaussian_clusters.png")
        plot_cfg = ClusterPlotConfig(; plot_kw..., boundaries = binfo, filename = png_path)
        t_plot = @elapsed fig = plot_clusters(smld, smld_out, plot_cfg)
        println("plot_clusters: ", round(t_plot, digits = 3), " s")

        # ── Cluster visualization (fixed-sigma "dot" render) ────────────────
        png_path2 = joinpath(results_dir, stem * "_dot_clusters.png")
        plot_cfg2 = ClusterPlotConfig(
            pixel_size_nm   = pixel_size_nm,
            colormap        = get(plot_kw, :colormap, :inferno),
            render_strategy = GaussianRender(use_localization_precision = false,
                                              fixed_sigma = 1.0),
            boundaries      = binfo,
            boundary_color  = :white,
            show_centroids  = false,
            filename        = png_path2,
        )
        t_plot2 = @elapsed fig2 = plot_clusters(smld, smld_out, plot_cfg2)
        println("plot_clusters (dot): ", round(t_plot2, digits = 3), " s")
    end
end

main()
