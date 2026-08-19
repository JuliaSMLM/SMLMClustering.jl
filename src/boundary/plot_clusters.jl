# `ClusterPlotConfig` is core (no plotting deps); `plot_clusters`'s real
# implementation lives in `SMLMClusteringFiguresExt` (loads when both
# `CairoMakie` and `SMLMRender` are present) — mirroring the edge-classify
# report/plot split in `src/viz.jl`.

## ─────────────────────────────────────────────────────────────────────────────
## Cluster-visualization subroutine
## ─────────────────────────────────────────────────────────────────────────────

"""
Configuration for `plot_clusters`.  All fields have sensible defaults so the
caller can override only what matters.
"""
Base.@kwdef struct ClusterPlotConfig
    # ── background render ──────────────────────────────────────────
    pixel_size_nm   :: Float64  = 10.0      # rendered pixel size (nm); smaller = finer detail
    colormap        :: Symbol   = :inferno  # colormap applied to the Gaussian density image
    scalebar        :: Bool     = false     # draw a physical scale bar on the rendered image
    scalebar_length :: Float64  = 1.0       # μm — explicit length
    render_strategy :: Any      = nothing   # background rendering strategy; nothing = SMLMRender's GaussianRender()
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

The rendered image uses SMLMRender's `GaussianRender` strategy (or
`cfg.render_strategy` if set) with the pixel size and colormap specified in
`cfg`.  Physical μm coordinates are recovered from the emitter data so that
boundary polygons and centroid markers align correctly with the rendered
image.

## Optional boundary overlay
If `cfg.boundaries` is a `ClusterBoundaryInfo` (returned by `boundary_clusters`),
each cluster's boundary polygon is drawn as a coloured closed curve on top of
the image.  Colours cycle through a small qualitative palette.

## Output
Returns the CairoMakie `Figure`.  If `cfg.filename !== nothing` the figure is
also saved to that path (format determined by file extension, e.g. ".png").

Requires the `SMLMClusteringFiguresExt` extension — load BOTH `CairoMakie`
and `SMLMRender` (`using CairoMakie, SMLMRender`).
"""
plot_clusters(args...; kwargs...) = error(
    "plot_clusters requires the SMLMClusteringFiguresExt extension — load BOTH " *
    "CairoMakie and SMLMRender (`using CairoMakie, SMLMRender`).")
