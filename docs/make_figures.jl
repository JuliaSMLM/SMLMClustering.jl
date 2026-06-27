# Standalone figure generator for the SMLMClustering docs.
#
# NOT wired into make.jl — run it manually under the heavy `examples/` project,
# which path-devs CairoMakie / SMLMData / SMLMClustering (so the docs build itself
# stays lightweight and fast):
#
#   julia --threads=auto --project=examples docs/make_figures.jl
#
# It writes PNGs into docs/src/assets/ (committed); the pages embed them. Synthetic
# point fields are hand-built from SMLMData emitter types (controllable Gaussian
# clusters + uniform noise) — clearer for illustrating clustering than the blinking
# simulator. Each figure is wrapped in try/catch and prints key numbers so one
# failure never aborts the whole run.

using CairoMakie
using SMLMData
using SMLMClustering
using Statistics
using Random

CairoMakie.activate!(type = "png")

const ASSETS = joinpath(@__DIR__, "src", "assets")
isdir(ASSETS) || mkpath(ASSETS)

const PALETTE = Makie.wong_colors()

# 0 = noise → translucent gray; 1..K → distinct opaque colors
function clustercolor(id::Integer)
    id <= 0 && return RGBAf(0.62, 0.62, 0.62, 0.35)
    c = PALETTE[mod1(id, length(PALETTE))]
    return RGBAf(c.r, c.g, c.b, 0.9)
end

xs(smld) = Float64[e.x for e in smld.emitters]
ys(smld) = Float64[e.y for e in smld.emitters]
ids(smld) = Int[e.id for e in smld.emitters]

emitter(x, y; id = 0) = SMLMData.Emitter2DFit{Float64}(
    x, y, 1000.0, 0.0, 0.01, 0.01, 0.0, 0.0;
    frame = 1, dataset = 1, track_id = 0, id = id)

function make_smld(emitters; fov = 5.0, px = 0.1)
    N = round(Int, fov / px)
    cam = SMLMData.IdealCamera(1:N, 1:N, px)
    return SMLMData.BasicSMLD(emitters, cam, 1, 1)
end

# Gaussian clusters (spec = (cx, cy, n, σ) in µm) + uniform noise.
function cluster_field(specs; fov = 5.0, n_noise = 500, seed = 1)
    Random.seed!(seed)
    em = SMLMData.Emitter2DFit{Float64}[]
    for (k, (cx, cy, n, σ)) in enumerate(specs)
        for _ in 1:n
            push!(em, emitter(cx + σ * randn(), cy + σ * randn(); id = k))
        end
    end
    for _ in 1:n_noise
        push!(em, emitter(rand() * fov, rand() * fov; id = 0))
    end
    return make_smld(em; fov = fov)
end

function blank_axis(gp; title = "")
    ax = Axis(gp; title = title, aspect = DataAspect())
    hidedecorations!(ax); hidespines!(ax)
    return ax
end

scatter_labels!(ax, x, y, labels; ms = 3) =
    scatter!(ax, x, y; color = [clustercolor(l) for l in labels], markersize = ms)

save_fig(name, fig) = (p = joinpath(ASSETS, name); save(p, fig; px_per_unit = 2); println("  ✓ wrote $p"); p)

# ---------------------------------------------------------------------------
# 1. Methods-overview comparison grid: one field, every labeling backend.
# ---------------------------------------------------------------------------
function fig_comparison_grid()
    println("[comparison_grid]")
    specs = [(1.2, 1.2, 220, 0.045), (3.7, 1.0, 160, 0.05), (2.5, 2.6, 300, 0.06),
             (1.0, 3.8, 130, 0.04), (3.9, 3.9, 200, 0.05), (2.6, 4.2, 90, 0.035)]
    smld = cluster_field(specs; fov = 5.0, n_noise = 500, seed = 7)
    x, y = xs(smld), ys(smld)

    # local-contrast seed/support for point-hysteresis
    _, lc = cluster_statistics(smld, LocalContrastFeature(density_k = 20, background_k = 200))
    ct = lc.extras[:contrast_per_emitter]; fn = lc.extras[:log_density_per_emitter]
    floor_f = quantile(filter(isfinite, fn), 0.3)
    seed = isfinite.(ct) .& isfinite.(fn) .& (ct .> 0.4) .& (fn .> floor_f)
    support = isfinite.(ct) .& (ct .> -0.1)

    runs = [
        ("DBSCAN",        () -> cluster(smld, DBSCANConfig(eps_nm = 60.0, min_points = 10))),
        ("HDBSCAN",       () -> cluster(smld, HDBSCANConfig(min_points = 8, min_cluster_size = 20))),
        ("Hierarchical",  () -> cluster(smld, HierarchicalConfig(n_clusters = length(specs), linkage = :ward))),
        ("Voronoi",       () -> cluster(smld, VoronoiConfig(density_factor = 2.0, min_points = 12))),
        ("MRF density",   () -> cluster(smld, MRFDensityClusterConfig(min_points = 12, density_estimator = :knn, density_k = 20))),
        ("Point hyst.",   () -> cluster(smld, PointHysteresisConfig(graph_k = 12, min_points = 20); seed = seed, support = support)),
    ]

    fig = Figure(size = (1080, 720))
    Label(fig[0, 1:3], "The same localization field, six backends"; fontsize = 20, font = :bold)
    for (i, (name, run)) in enumerate(runs)
        r, c = fldmod1(i, 3)
        ax = blank_axis(fig[r, c])
        try
            out, info = run()
            scatter_labels!(ax, x, y, ids(out))
            ax.title = "$name — $(info.n_clusters) clusters"
            println("  $name: n_clusters=$(info.n_clusters), n_noise=$(info.n_noise)")
        catch err
            scatter_labels!(ax, x, y, zeros(Int, length(x)))
            ax.title = "$name — failed"
            println("  $name FAILED: $err")
        end
    end
    save_fig("comparison_grid.png", fig)
end

# ---------------------------------------------------------------------------
# 2. MRF density-regime: continuous per-emitter density → 3 spatially-coherent regimes.
# ---------------------------------------------------------------------------
function fig_mrf_regimes()
    println("[mrf_regimes]")
    # Three genuine density regimes: sparse background, a medium-density extended
    # patch, and dense cores embedded in it — the multi-regime case a single ε can't fit.
    Random.seed!(21)
    fov = 5.0
    em = SMLMData.Emitter2DFit{Float64}[]
    for _ in 1:500                                            # sparse background
        push!(em, emitter(rand() * fov, rand() * fov; id = 0))
    end
    for _ in 1:1000                                           # medium extended patch
        push!(em, emitter(2.5 + 0.6 * randn(), 2.5 + 0.6 * randn(); id = 0))
    end
    for (cx, cy) in [(1.9, 2.4), (2.6, 3.1), (3.1, 2.2)]      # dense cores within it
        for _ in 1:220
            push!(em, emitter(cx + 0.04 * randn(), cy + 0.04 * randn(); id = 0))
        end
    end
    smld = make_smld(em; fov = fov)
    x, y = xs(smld), ys(smld)

    _, vd = cluster_statistics(smld, VoronoiDensityConfig())
    dens = vd.extras[:density_per_emitter]
    logd = [d > 0 ? log10(d) : NaN for d in dens]

    # use the default :voronoi estimator so the regimes are binned from exactly the
    # per-emitter density shown in the left panel
    out, info = cluster(smld, MRFDensityClusterConfig(
        n_regimes = 3, density_estimator = :voronoi, min_points = 15))
    reg = out.metadata["mrf_regime_per_emitter"]

    fig = Figure(size = (900, 460))
    Label(fig[0, 1:2], "MRF density-regime: continuous density → coherent regimes";
          fontsize = 18, font = :bold)
    ax1 = blank_axis(fig[1, 1]; title = "per-emitter density (log₁₀ µm⁻²)")
    finite = isfinite.(logd)
    scatter!(ax1, x[finite], y[finite]; color = logd[finite], colormap = :viridis, markersize = 3)
    ax2 = blank_axis(fig[1, 2]; title = "MRF regimes (background / patch / dense cores)")
    scatter_labels!(ax2, x, y, reg)
    for r in 0:3
        println("  regime $r: $(count(==(r), reg)) emitters")
    end
    println("  clusters (foreground CC): $(info.n_clusters)")
    save_fig("mrf_regimes.png", fig)
end

# ---------------------------------------------------------------------------
# 3. Hopkins: clustered (H→1) vs uniform random (H≈0.5).
# ---------------------------------------------------------------------------
function fig_hopkins()
    println("[hopkins]")
    clustered = cluster_field([(1.3, 1.3, 300, 0.06), (3.6, 1.4, 300, 0.06),
                               (2.4, 3.6, 300, 0.06)]; n_noise = 150, seed = 3)
    Random.seed!(5)
    uniform = make_smld([emitter(rand() * 5, rand() * 5; id = 0) for _ in 1:1050])

    Hc = cluster_statistics(clustered, HopkinsConfig(n_samples = 100, random_repeats = 20, seed = 1))[2].statistic
    Hu = cluster_statistics(uniform,   HopkinsConfig(n_samples = 100, random_repeats = 20, seed = 1))[2].statistic
    println("  H_clustered=$(round(Hc, digits=3))  H_uniform=$(round(Hu, digits=3))")

    fig = Figure(size = (760, 420))
    Label(fig[0, 1:2], "Hopkins clustering tendency"; fontsize = 18, font = :bold)
    ax1 = blank_axis(fig[1, 1]; title = "clustered   H = $(round(Hc, digits=2))")
    scatter!(ax1, xs(clustered), ys(clustered); color = (:firebrick, 0.6), markersize = 3)
    ax2 = blank_axis(fig[1, 2]; title = "uniform random   H = $(round(Hu, digits=2))")
    scatter!(ax2, xs(uniform), ys(uniform); color = (:steelblue, 0.6), markersize = 3)
    save_fig("hopkins.png", fig)
end

# ---------------------------------------------------------------------------
# 4. Voronoi density: per-emitter density feature → SR-Tesseler clusters.
# ---------------------------------------------------------------------------
function fig_voronoi()
    println("[voronoi]")
    smld = cluster_field([(1.3, 1.3, 320, 0.05), (3.5, 1.5, 240, 0.05),
                          (2.5, 3.5, 360, 0.06), (3.9, 3.8, 180, 0.045)];
                         n_noise = 500, seed = 9)
    x, y = xs(smld), ys(smld)
    _, vd = cluster_statistics(smld, VoronoiDensityConfig())
    dens = vd.extras[:density_per_emitter]
    logd = [d > 0 ? log10(d) : NaN for d in dens]
    out, info = cluster(smld, VoronoiConfig(density_factor = 2.0, min_points = 12))

    fig = Figure(size = (820, 420))
    Label(fig[0, 1:2], "Voronoi / SR-Tesseler"; fontsize = 18, font = :bold)
    ax1 = blank_axis(fig[1, 1]; title = "per-emitter Voronoi density (log₁₀ µm⁻²)")
    finite = isfinite.(logd)
    scatter!(ax1, x[finite], y[finite]; color = logd[finite], colormap = :magma, markersize = 3)
    ax2 = blank_axis(fig[1, 2]; title = "dense-cell clusters — $(info.n_clusters) found")
    scatter_labels!(ax2, x, y, ids(out))
    save_fig("voronoi_density.png", fig)
end

# ---------------------------------------------------------------------------
# 5. Local contrast: cancels a density gradient that raw density cannot.
# ---------------------------------------------------------------------------
function fig_local_contrast()
    println("[local_contrast]")
    Random.seed!(13)
    fov = 5.0
    em = SMLMData.Emitter2DFit{Float64}[]
    # background whose density rises left→right (rejection sampling on x)
    nbg = 6000
    while count(e -> e.id == 0, em) < nbg
        x = rand() * fov
        rand() < (0.15 + 0.85 * x / fov) || continue
        push!(em, emitter(x, rand() * fov; id = 0))
    end
    # a few compact clusters across the gradient
    for (k, (cx, cy)) in enumerate([(1.0, 2.5), (2.5, 2.5), (4.0, 2.5)])
        for _ in 1:220
            push!(em, emitter(cx + 0.05 * randn(), cy + 0.05 * randn(); id = k))
        end
    end
    smld = make_smld(em; fov = fov)
    x, y = xs(smld), ys(smld)
    _, lc = cluster_statistics(smld, LocalContrastFeature(density_k = 25, background_k = 400))
    fine = lc.extras[:log_density_per_emitter]
    cont = lc.extras[:contrast_per_emitter]

    fig = Figure(size = (1180, 380))
    Label(fig[0, 1:2], "Local contrast cancels a baseline density gradient"; fontsize = 18, font = :bold)
    ax1 = blank_axis(fig[1, 1]; title = "raw kNN log-density — gradient dominates")
    f1 = isfinite.(fine)
    scatter!(ax1, x[f1], y[f1]; color = fine[f1], colormap = :viridis, markersize = 3)
    ax2 = blank_axis(fig[1, 2]; title = "local contrast — clusters pop out")
    f2 = isfinite.(cont)
    scatter!(ax2, x[f2], y[f2]; color = cont[f2], colormap = :viridis,
             colorrange = (quantile(cont[f2], 0.02), quantile(cont[f2], 0.98)), markersize = 3)
    save_fig("local_contrast.png", fig)
end

# ---------------------------------------------------------------------------
# 6. Edge classification: outside / membrane / interior on a synthetic cell.
# ---------------------------------------------------------------------------
function fig_edge_classify()
    println("[edge_classify]")
    Random.seed!(17)
    fov = 5.0; cx, cy, R = 2.5, 2.5, 1.5
    em = SMLMData.Emitter2DFit{Float64}[]
    for _ in 1:1700                                   # interior fill
        r = R * 0.92 * sqrt(rand()); θ = 2π * rand()
        push!(em, emitter(cx + r * cos(θ), cy + r * sin(θ)))
    end
    for _ in 1:1100                                   # denser membrane ring
        θ = 2π * rand(); rr = R + 0.03 * randn()
        push!(em, emitter(cx + rr * cos(θ), cy + rr * sin(θ)))
    end
    for _ in 1:500                                    # sparse uniform field background
        push!(em, emitter(rand() * fov, rand() * fov))
    end
    smld = make_smld(em; fov = fov)
    _, info = classify_emitters(smld, OuterPolygonConfig())
    x, y = xs(smld), ys(smld)
    colmap = Dict(:outside => RGBAf(0.7, 0.7, 0.7, 0.3),
                  :membrane => RGBAf(0.90, 0.49, 0.13, 0.9),
                  :interior => RGBAf(0.20, 0.49, 0.72, 0.7))
    fig = Figure(size = (560, 540))
    ax = blank_axis(fig[1, 1]; title = "Edge classification (OuterPolygon)")
    scatter!(ax, x, y; color = [colmap[c] for c in info.class], markersize = 3)
    poly = info.outer_polygon                       # draw the alpha-shape outer loop
    if !isempty(poly)
        plx = Float64[p[1] for p in poly]; ply = Float64[p[2] for p in poly]
        push!(plx, plx[1]); push!(ply, ply[1])      # close the loop
        lines!(ax, plx, ply; color = :black, linewidth = 2)
    end
    n_o = count(==(:outside), info.class); n_m = count(==(:membrane), info.class); n_i = count(==(:interior), info.class)
    println("  outside=$n_o membrane=$n_m interior=$n_i")
    save_fig("edge_classify.png", fig)
end

function main()
    println("Generating SMLMClustering doc figures → $ASSETS")
    for f in (fig_comparison_grid, fig_mrf_regimes, fig_hopkins,
              fig_voronoi, fig_local_contrast, fig_edge_classify)
        try
            f()
        catch err
            println("  !! $(nameof(f)) FAILED: $err")
        end
    end
    println("done.")
end

main()
