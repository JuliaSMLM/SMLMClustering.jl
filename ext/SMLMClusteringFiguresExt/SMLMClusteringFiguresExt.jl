module SMLMClusteringFiguresExt

using SMLMClustering
using SMLMData
using CairoMakie
using SMLMRender

# ---- class render (SMLMRender) ----------------------------------------------

# class -> single render color (genmab-canonical)
const _CLASS_RENDER = ((:outside, :gray70), (:membrane, :gold), (:interior, :steelblue))

"""
    render_classes(smld, class; extra=[], zoom=20, radius_factor=20, line_width=1,
                   filename=nothing)

SMLMRender circle render of `smld` split by `class` (:outside gray / :membrane gold /
:interior steelblue): each class is a separately-rendered layer (so each
self-normalizes) composed additively onto one shared target. `extra` adds layers
`(mask, color, name)` — e.g. `[(structure_id .!= 0, :red, "structure")]` for a
4-class image; `extra=[]` gives the 3-class default. Saves to `filename` if given;
returns the composed image.

Requires both `CairoMakie` and `SMLMRender` loaded (this `SMLMClusteringFiguresExt`).
"""
function SMLMClustering.render_classes(smld, class::AbstractVector;
        extra = [], zoom = 20, radius_factor = 20, line_width = 1, filename = nothing)
    length(class) == length(smld.emitters) ||
        throw(ArgumentError("class length $(length(class)) ≠ n emitters $(length(smld.emitters))"))
    target = create_target_from_smld(smld; zoom = zoom)
    strat = CircleRender(; radius_factor = radius_factor, line_width = line_width,
                         use_localization_precision = true)
    layers = Tuple{Any, Any}[]
    for (cl, col) in _CLASS_RENDER
        push!(layers, (class .== cl, col))
    end
    for (mask, col, _name) in extra
        push!(layers, (collect(Bool, mask), col))
    end
    imgs = Any[]
    for (mask, col) in layers
        any(mask) || continue
        sub = SMLMData.BasicSMLD(smld.emitters[mask], smld.camera,
                                 smld.n_frames, smld.n_datasets)
        (img, _) = render(sub; strategy = strat, color = col, target = target,
                          clip_percentile = nothing)
        push!(imgs, img)
    end
    isempty(imgs) && throw(ArgumentError("render_classes: no emitters to render"))
    composed = length(imgs) == 1 ? imgs[1] : compose(imgs...)
    filename === nothing || save_image(filename, composed)
    return composed
end

# ---- CairoMakie overlay + class-fraction panel ------------------------------

const _CLASS_STYLE = ((:interior, (:steelblue, 0.30)),
                      (:membrane, (:orange, 0.85)),
                      (:outside,  (:firebrick, 0.40)))

function _ring!(ax, ring; kwargs...)
    isempty(ring) && return
    px = Float64[p[1] for p in ring]; py = Float64[p[2] for p in ring]
    push!(px, px[1]); push!(py, py[1])
    lines!(ax, px, py; kwargs...)
end

function _overlay(report, path)
    info = report.info
    fig = Figure(size = (900, 900))
    ax = Axis(fig[1, 1]; title = "edge mask — $(method_name(info.config))", aspect = DataAspect())
    hidedecorations!(ax); hidespines!(ax)
    for (cl, col) in _CLASS_STYLE
        idx = findall(==(cl), info.class)
        isempty(idx) || scatter!(ax, report.x_um[idx], report.y_um[idx];
                                 color = col, markersize = 2.0, label = String(cl))
    end
    for cell in info.cells
        _ring!(ax, cell.outer; color = :black, linewidth = 2.0)
        for h in cell.holes
            _ring!(ax, h; color = (:black, 0.7), linewidth = 1.2, linestyle = :dash)
        end
    end
    axislegend(ax; position = :rt, framevisible = false)
    save(path, fig)
end

function _fractions(report, path)
    fr = report.fractions
    fig = Figure(size = (520, 400))
    ax = Axis(fig[1, 1]; title = "class fractions", ylabel = "fraction",
              xticks = (1:3, ["interior", "membrane", "outside"]))
    barplot!(ax, 1:3, [fr.interior, fr.membrane, fr.outside];
             color = [:steelblue, :orange, :firebrick])
    ylims!(ax, 0, 1)
    save(path, fig)
end

"""
    plot_edge_report(report; output_dir, zoom_overlay=10, zoom_render=20, prefix="edge")
        -> Vector{String}

Write the standard edge-mask figure series into `output_dir`, returning the saved
paths: `<prefix>_classes.png` (SMLMRender class image), `<prefix>_overlay.png`
(CairoMakie polygon overlay over class-colored localizations) and
`<prefix>_fractions.png` (class-fraction bar). `zoom_render` sets the class-render
resolution; `zoom_overlay` is reserved for a future render-backed overlay variant.

Requires both `CairoMakie` and `SMLMRender` loaded (this `SMLMClusteringFiguresExt`).
"""
function SMLMClustering.plot_edge_report(report::SMLMClustering.EdgeReport;
        output_dir::AbstractString, zoom_overlay = 10, zoom_render = 20, prefix = "edge")
    mkpath(output_dir)
    paths = String[]
    pr = joinpath(output_dir, "$(prefix)_classes.png")
    SMLMClustering.render_classes(report.smld, report.info.class; zoom = zoom_render, filename = pr)
    push!(paths, pr)
    po = joinpath(output_dir, "$(prefix)_overlay.png"); _overlay(report, po); push!(paths, po)
    pf = joinpath(output_dir, "$(prefix)_fractions.png"); _fractions(report, pf); push!(paths, pf)
    return paths
end

end # module
