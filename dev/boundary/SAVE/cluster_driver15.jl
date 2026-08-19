using SMLMAnalysis        # render, load_smld, RenderConfig, GaussianRender, BasicSMLD, …
using SMLMData            # filter_roi
using SMLMClustering      # cluster, cluster_statistics, boundary_clusters, …
using CairoMakie          # Figure, Axis, image!, lines!, scatter!, save
using Statistics          # mean, std

## ─────────────────────────────────────────────────────────────────────────────
## Cluster-visualization subroutine
## ─────────────────────────────────────────────────────────────────────────────

"""
Configuration for `plot_clusters`.  All fields have sensible defaults so the
caller can override only what matters.
"""
Base.@kwdef struct ClusterPlotConfig
    # ── background render ──────────────────────────────────────────
    pixel_size_nm   :: Float64  = 100.0     # rendered pixel size (nm); smaller = finer detail
    colormap        :: Symbol   = :inferno  # colormap applied to the Gaussian density image
#   scalebar        :: Bool     = true      # draw a physical scale bar on the rendered image
    scalebar        :: Bool     = false     # draw a physical scale bar on the rendered image
    scalebar_length :: Float64  = 1.0       # μm — explicit length
    render_strategy :: Any      = GaussianRender()  # background rendering strategy
    color_by        :: Union{Symbol, Nothing} = nothing  # field to color localizations by (e.g. :id)
    categorical     :: Bool     = false     # treat color_by field as categorical

    # ── cluster centroid markers ───────────────────────────────────
    show_centroids :: Bool    = true    # draw ×  at each cluster's centroid
    centroid_size  :: Float64 = 12.0    # marker size in screen points
    centroid_color :: Any     = :white  # marker colour (any Makie-compatible colour)

    # ── boundary-polygon overlay (optional) ───────────────────────
    # Pass the ClusterBoundaryInfo returned by boundary_clusters() to draw
    # each cluster's alpha-shape outline.  nothing = no boundary overlay.
    boundaries    :: Union{ClusterBoundaryInfo, Nothing} = nothing
    boundary_linewidth :: Float64 = 1.0  # outline line width in points
    boundary_color :: Union{Symbol, Nothing} = nothing  # nothing = cycle palette; else fixed color for all boundaries

    # ── output ────────────────────────────────────────────────────
    filename    :: Union{String, Nothing} = nothing  # PNG save path; nothing = skip saving
    figure_size :: Tuple{Int,Int}         = (900, 900)  # figure size in pixels
end

"""
    plot_clusters(smld, smld_out, cfg) -> Figure

Render a Gaussian super-resolution image of `smld` (all localizations) and
overlay cluster indicators derived from `smld_out` (emitter `.id` field).

The rendered image uses SMLMRender's `GaussianRender` strategy with the pixel
size and colormap specified in `cfg`.  Physical μm coordinates are recovered
from the emitter data so that boundary polygons and centroid markers align
correctly with the rendered image.

## Optional boundary overlay
If `cfg.boundaries` is a `ClusterBoundaryInfo` (returned by `boundary_clusters`),
each cluster's boundary polygon is drawn as a coloured closed curve on top of
the image.  Colours cycle through a small qualitative palette.

## Output
Returns the CairoMakie `Figure`.  If `cfg.filename !== nothing` the figure is
also saved to that path (format determined by file extension, e.g. ".png").
"""
function plot_clusters(smld     :: BasicSMLD,
                       smld_out :: BasicSMLD,
                       cfg      :: ClusterPlotConfig)

    # ── 1. Render the Gaussian background image ───────────────────────────────
    # render() in pixel_size mode derives data bounds from emitter positions
    # and adds a 5 % margin on each side (SMLMRender default).
    (bg_img, rinfo) = render(smld, RenderConfig(
        strategy        = cfg.render_strategy,
#       pixel_size      = cfg.pixel_size_nm,
        pixel_size      = cfg.pixel_size_nm / 10,
        colormap        = cfg.colormap,
        scalebar        = cfg.scalebar,
        scalebar_length = cfg.scalebar_length,
        color_by        = cfg.color_by,
        categorical     = cfg.categorical,
    ))

    # ── 2. Re-derive the physical coordinate bounds ───────────────────────────
    # We mimic SMLMRender's create_target_from_smld(; pixel_size) so that the
    # CairoMakie axis uses the same μm ranges as the rendered pixels.
    emitter_x = [e.x for e in smld.emitters]
    emitter_y = [e.y for e in smld.emitters]
    x_min_d, x_max_d = extrema(emitter_x)
    y_min_d, y_max_d = extrema(emitter_y)
    mg = 0.05                                  # 5 % margin (matches SMLMRender default)
    x_min = x_min_d - mg * (x_max_d - x_min_d)
    x_max = x_max_d + mg * (x_max_d - x_min_d)
    y_min = y_min_d - mg * (y_max_d - y_min_d)
    y_max = y_max_d + mg * (y_max_d - y_min_d)

    # ── 3. Assemble the figure ────────────────────────────────────────────────
    n_clustered = count(e -> e.id > 0, smld_out.emitters)
    n_total     = length(smld_out.emitters)

    fig = Figure(size = cfg.figure_size)
    ax  = Axis(fig[1, 1];
        xlabel    = "x (μm)",
        ylabel    = "y (μm)",
        title     = "Clusters: $n_clustered / $n_total localizations clustered",
        yreversed = true,   # SMLM y increases downward (camera convention)
    )

    # Display the rendered image.
    # bg_img has layout [row, col] = [y_pixel, x_pixel].
    # CairoMakie's image! expects data[i,j] at (x[i], y[j]), so we transpose
    # to [col, row] = [x_pixel, y_pixel] before passing.
    image!(ax, (x_min, x_max), (y_min, y_max),
           permutedims(bg_img, (2, 1)))

    # ── 4. Cluster centroid markers ───────────────────────────────────────────
    if cfg.show_centroids
        # Group clustered emitters by (dataset, cluster_id) — id == 0 is noise.
        clustered = filter(e -> e.id > 0, smld_out.emitters)
        if !isempty(clustered)
            cluster_keys = sort!(unique((e.dataset, e.id) for e in clustered))
            for (ds, cid) in cluster_keys
                ex = [e.x for e in clustered if e.dataset == ds && e.id == cid]
                ey = [e.y for e in clustered if e.dataset == ds && e.id == cid]
                # Centroid = arithmetic mean of cluster member positions.
                cx = sum(ex) / length(ex)
                cy = sum(ey) / length(ey)
                scatter!(ax, [cx], [cy];
                    marker     = :xcross,
                    markersize = cfg.centroid_size,
                    color      = cfg.centroid_color,
                )
            end
        end
    end

    # ── 5. Boundary-polygon overlay (optional) ────────────────────────────────
    if cfg.boundaries !== nothing
        binfo = cfg.boundaries

        # Qualitative colour palette — cycles for data sets with many clusters.
        palette = [:cyan, :yellow, :magenta, :lime, :orange,
                   :red,  :blue,   :white,   :pink, :aqua]

        for (j, _key) in enumerate(binfo.cluster_keys)
            bx = binfo.boundaries_x[j]   # closed polygon x-coords in μm
            by = binfo.boundaries_y[j]   # closed polygon y-coords in μm
            isempty(bx) && continue       # skip clusters that produced no boundary

            line_color = cfg.boundary_color !== nothing ? cfg.boundary_color :
                         palette[mod1(j, length(palette))]
            lines!(ax, bx, by;
                color     = line_color,
                linewidth = cfg.boundary_linewidth,
            )
        end
    end

    # ── 6. Save to file if a path was provided ────────────────────────────────
    if cfg.filename !== nothing
        save(cfg.filename, fig)
        println("Saved → ", cfg.filename)
    end

    return fig
end

## ─────────────────────────────────────────────────────────────────────────────
## Main driver
## ─────────────────────────────────────────────────────────────────────────────

# Path to input h5 smld dataset
#smld_path = "/mnt/nas/cellpath/papers-geometry/julia-analysis-results/dstorm/20231218_IgE-DF3/DF310nM5min/Cell_06/14_render/gaussianrender_tab10_id_50x.png"
smld_path = "/mnt/nas/lidkelab/Personal Folders/MJW/Julia/publish/resultsH5/" *
            "Data_2025-11-20-10-52-19_result_smld.h5"
smld = load_smld(smld_path)
smld = filter_roi(smld, (12.5, 17.5), (12.5, 17.5))   # 5 μm × 5 μm centered at (15, 15) μm
#smld = filter_roi(smld, (12.5, 13.5), (12.5, 13.5))

# Path to directory to place results.
results_dir = "/mnt/nas/lidkelab/Personal Folders/MJW/Julia/"

# Cluster within each dataset independently; must match in the configuration
# parameters for cluster and boundary_clusters,
per_dataset = false

# ── Clustering (DBSCAN) ───────────────────────────────────────────────────────
cfg_dbscan = DBSCANConfig(
    eps_nm             = 25.0,  # neighborhood radius in nm (required)
    min_points         = 5,     # core-point threshold / min cluster size
    use_3d             = false, # include z-coordinate
    per_dataset        = per_dataset, # cluster within each dataset independently
    remove_unclustered = false,
)
(smld_out, info_dbscan) = cluster(smld, cfg_dbscan)
println(smld_out)
println(info_dbscan)

# ── Cluster boundaries (alpha-shape) ─────────────────────────────────────────
cfg_bnd = BoundaryClustersConfig(
#   shrink      = 0.5,  # 0 = convex hull, 1 = tightest alpha-shape
    shrink      = 0.05,  # 0 = convex hull, 1 = tightest alpha-shape
    per_dataset = per_dataset, # must match the per_dataset value used in cluster()
)
(_, binfo) = boundary_clusters(smld_out, cfg_bnd)
println(binfo)
# binfo.boundaries_x[j] / binfo.boundaries_y[j] trace the boundary polygon of
# cluster binfo.cluster_keys[j] in μm (closed: first point == last point).

# ── Cluster boundary statistics ───────────────────────────────────────────────
(_, sinfo) = boundary_statistics(smld_out, binfo)
println(sinfo)

valid_areas  = filter(!isnan, sinfo.areas)
valid_perims = filter(!isnan, sinfo.perimeters)

if !isempty(valid_areas)
    n_a = length(valid_areas)
    println("Area:      = ", round(mean(valid_areas), digits = 4),
            "  +/- ", round(std(valid_areas) / sqrt(n_a), digits = 4), " μm²",
            "  (n = ", n_a, ")")
end
if !isempty(valid_perims)
    n_p = length(valid_perims)
    println("Perimeter: = ", round(mean(valid_perims), digits = 4),
            "  +/- ", round(std(valid_perims) / sqrt(n_p), digits = 4), " μm ",
            "  (n = ", n_p, ")")
end

# ── Cluster visualization ─────────────────────────────────────────────────────
# Derive a PNG filename from the data file (placed next to this script).
png_path = results_dir * splitext(basename(smld_path))[1] * "_Gaussian_clusters.png"
plot_cfg = ClusterPlotConfig(
    pixel_size_nm = 97.8 / 10,     # 10 nm pixels → fine enough for 50 nm DBSCAN
    colormap      = :inferno,      # good contrast for SMLM density images
    color_by      = :id,           # color each localization by its cluster id
    categorical   = true,
    boundaries    = binfo,         # draw alpha-shape outlines for each cluster
    show_centroids = false,
    filename      = png_path,
)
t_plot = @elapsed fig = plot_clusters(smld, smld_out, plot_cfg)
println("plot_clusters: ", round(t_plot, digits = 3), " s")

# ── Cluster visualization (fixed-sigma "dot" render) ──────────────────────────
png_path2 = results_dir * splitext(basename(smld_path))[1] * "_dot_clusters.png"
plot_cfg2 = ClusterPlotConfig(
    pixel_size_nm   = 97.8 / 10,
    colormap        = :inferno,
    render_strategy = GaussianRender(use_localization_precision = false,
                                      fixed_sigma = 1.0),
    boundaries      = binfo,
    boundary_color  = :white,
    show_centroids  = false,
    filename        = png_path2,
)
t_plot2 = @elapsed fig2 = plot_clusters(smld, smld_out, plot_cfg2)
println("plot_clusters (dot): ", round(t_plot2, digits = 3), " s")
