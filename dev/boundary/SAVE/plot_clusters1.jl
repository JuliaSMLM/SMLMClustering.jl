using SMLMAnalysis        # render, RenderConfig, GaussianRender, BasicSMLD
using SMLMClustering      # ClusterBoundaryInfo
using CairoMakie          # Figure, Axis, image!, lines!, scatter!, save

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
