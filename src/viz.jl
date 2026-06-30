"""
Edge-mask reporting + figures. The *compute* half is core (no plotting deps); the
*plot* half lives in `SMLMClusteringFiguresExt` (loads when both `CairoMakie` and
`SMLMRender` are present) — mirroring SMLMBaGoL's report/plot split.

- [`compute_edge_report`](@ref)`(smld, info)` -> [`EdgeReport`](@ref)   (core)
- [`write_edge_report`](@ref)`(report; output_dir)`                      (core; text/TSV diagnostics)
- `plot_edge_report(report; output_dir, zoom_overlay, zoom_render, prefix)` -> `Vector{String}`
  and `render_classes(smld, class; extra, zoom)`                        (extension)
- [`class_codes`](@ref)`(info)` -> `Vector{Int}`                         (core)
"""

"""
    EdgeReport

Figure-data derivative of an [`EdgeClassifyInfo`](@ref): the classified SMLD plus
per-emitter coordinates/classes and per-class fractions needed to write diagnostics
and render the standard edge figures. Produced by [`compute_edge_report`](@ref);
consumed by [`write_edge_report`](@ref) and `plot_edge_report`.
"""
struct EdgeReport{S, C<:AbstractEdgeClassifyConfig}
    smld::S
    info::EdgeClassifyInfo{C}
    x_um::Vector{Float64}
    y_um::Vector{Float64}
    fractions::NamedTuple{(:interior, :membrane, :outside), NTuple{3, Float64}}
end

"""
    compute_edge_report(smld, info::EdgeClassifyInfo) -> EdgeReport

Bundle the classified SMLD and `info` into the figure-data report. Pass the SMLD
returned by `classify_emitters` — its emitters are unchanged, so `info.class` is
1:1 with `smld.emitters`.
"""
function compute_edge_report(smld, info::EdgeClassifyInfo)
    n = length(smld.emitters)
    n == length(info.class) ||
        throw(ArgumentError("smld has $n emitters but info classifies $(length(info.class))"))
    x = Vector{Float64}(undef, n); y = Vector{Float64}(undef, n)
    @inbounds for i in 1:n
        x[i] = smld.emitters[i].x; y[i] = smld.emitters[i].y
    end
    tot = max(n, 1)
    fr = (interior = info.n_interior / tot,
          membrane = info.n_membrane / tot,
          outside  = info.n_outside  / tot)
    return EdgeReport(smld, info, x, y, fr)
end

"""
    write_edge_report(report; output_dir, condition="", cell="") -> String

Write the edge-classify text/TSV diagnostics (classified.tsv, polygon_loops.tsv,
loop_diagnostics.csv, params.json, manifest.json) for `report` into `output_dir`
(folds [`write_edge_artifacts`](@ref)). Returns `output_dir`.
"""
function write_edge_report(report::EdgeReport; output_dir::AbstractString,
                           condition::AbstractString = "", cell::AbstractString = "")
    meta = hasproperty(report.smld, :metadata) ? report.smld.metadata : nothing
    if meta isa AbstractDict
        # Drop the classify OUTPUTS that classify_emitters wrote into the SMLD
        # (edge_cells is a Vector{CellPolygon}, not JSON-serializable provenance);
        # keep only the genuine input metadata for params.json.
        meta = Dict{String,Any}(string(k) => v for (k, v) in meta if !startswith(string(k), "edge_"))
    end
    write_edge_artifacts(output_dir, report.info, report.x_um, report.y_um;
                         condition = condition, cell = cell, smld_input_meta = meta)
    return output_dir
end

# Extension entry points. The real methods live in `SMLMClusteringFiguresExt`, which
# loads only when BOTH `CairoMakie` and `SMLMRender` are present. These fallbacks give a
# clear message instead of a bare `MethodError` when the extension isn't active.
const _FIG_EXT_HINT = "requires the SMLMClusteringFiguresExt extension — load BOTH " *
                      "CairoMakie and SMLMRender (`using CairoMakie, SMLMRender`)."
plot_edge_report(args...; kwargs...) = error("plot_edge_report " * _FIG_EXT_HINT)
render_classes(args...; kwargs...)   = error("render_classes "   * _FIG_EXT_HINT)

"""
    class_codes(info::EdgeClassifyInfo) -> Vector{Int}

Per-emitter integer class code for categorical rendering: `outside = 0`,
`membrane = 1`, `interior = 2`. The `0` for outside matches SMLMRender's
reserved-gray id.
"""
function class_codes(info::EdgeClassifyInfo)
    codes = Vector{Int}(undef, length(info.class))
    @inbounds for i in eachindex(info.class)
        c = info.class[i]
        codes[i] = c === :interior ? 2 : (c === :membrane ? 1 : 0)
    end
    return codes
end
