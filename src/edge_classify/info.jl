"""
Result records for `classify_emitters`.

`EdgeClassifyInfo <: SMLMData.AbstractSMLMInfo` is the sibling of `ClusterInfo` /
`ClusterStatisticsInfo`. Because emitters have no categorical class field (unlike
cluster `id`s), the canonical per-emitter answer lives here (`info.class`, with the
`in_cell` / `interior_mask` accessors). It is deliberately *not* mirrored into the
SMLD metadata — a per-emitter side-list would desync the moment a downstream step
subsets emitters; consume it at the classify point. Only the per-cell mask GEOMETRY
travels in metadata (`edge_cells`, `edge_outer_polygon`).
"""

# Valid per-emitter class symbols (serialized as Strings on disk).
const EDGE_CLASSES = (:outside, :membrane, :interior)

"""
    LoopDiagnostic

Per-loop diagnostic record (one row of `loop_diagnostics.csv`).
"""
struct LoopDiagnostic
    loop_id::Int
    vertex_count::Int
    area_um2::Float64
    n_emitters_inside::Int
    frac_in_fov::Float64
    frac_dense::Float64
    median_rhoK::Float64
    used_in_outer::Bool
    heuristic_type::String
end

"""
    ConcavityMetricReport

Boundary-proximal concavity-error metric for the outer-polygon classifier (see
[`compute_concavity_metric`](@ref)).
"""
struct ConcavityMetricReport
    buffer_um::Float64
    asym_R_nm::Float64
    asym_gate::Float64
    rho_lo::Float64
    n_interior::Int
    n_eligible::Int
    n_suspect::Int
    n_suspect_interior_fov::Int
    n_suspect_fov_edge::Int
    suspect_x_um::Vector{Float64}
    suspect_y_um::Vector{Float64}
    suspect_is_fov_edge::BitVector
end

"""
    EdgeClassifyInfo{C} <: SMLMData.AbstractSMLMInfo

Result of `classify_emitters`. `class` is the authoritative per-emitter answer
(`:outside`, `:membrane`, `:interior`), a partition of the input set.

`cells::MultiCellMask` is the published mask: one `CellPolygon` (an outer ring plus
optional internal holes) per distinct cell in the FOV, ordered largest-first. It is
the drawn boundary and the Hopkins `region = :metadata` observation window.
`outer_polygon` is `cells[1].outer` (the dominant cell's outer ring), retained for
back-compat.

`inside_outer` is geometric membership in the mask (`in_region(cells)`), and
`dist_to_outer_um` is the distance to the nearest **real** boundary segment
(`NaN` when not inside). Boundary segments lying on a truncated FOV edge are
excluded from that distance, so a field-of-view cut is never labeled membrane —
`membrane` is the band within `membrane_nm` of a real cell edge. `in_cell(info)`
(`class .!= :outside`) equals `inside_outer`. Read `class` for the answer.

`loops` are the raw alpha-shape boundary loops (serialized to `polygon_loops.tsv`);
`config` holds the concrete config that ran (honest provenance).
"""
struct EdgeClassifyInfo{C<:AbstractEdgeClassifyConfig} <: SMLMData.AbstractSMLMInfo
    n_emitters::Int
    class::Vector{Symbol}
    inside_outer::BitVector
    dist_to_outer_um::Vector{Float64}
    outer_polygon::Vector{NTuple{2,Float64}}
    cells::MultiCellMask
    loops::Vector{Vector{NTuple{2,Float64}}}
    loop_diagnostics::Vector{LoopDiagnostic}
    config::C
    fov_um::NTuple{4,Float64}
    truncated_sides::NamedTuple{(:L, :R, :B, :T), NTuple{4,Bool}}
    runtime_s::Float64
    n_outside::Int
    n_membrane::Int
    n_interior::Int
end

"""
    in_cell(info) -> BitVector

Topological cell membership, `== (info.class .!= :outside)`. Equals `inside_outer`
(both are mask membership via `in_region`).
"""
in_cell(info::EdgeClassifyInfo) = info.class .!= :outside

"""
    interior_mask(info) -> BitVector

Per-emitter interior mask, `== (info.class .== :interior)` — the strictly-interior
emitters (membrane and outside both excluded). Unlike [`in_cell`](@ref) (which keeps
membrane), this is the subset a downstream analysis usually carries forward, and is
the boolean to AND with any other per-emitter mask (e.g. a separate structure mask)
before subsetting. `classify_emitters` stays non-destructive — it never drops
emitters — so the caller composes masks and subsets once.
"""
interior_mask(info::EdgeClassifyInfo) = info.class .== :interior

"""
    interior_fraction(info) -> Float64
"""
interior_fraction(info::EdgeClassifyInfo) = info.n_interior / max(info.n_emitters, 1)

function Base.show(io::IO, info::EdgeClassifyInfo)
    print(io, "EdgeClassifyInfo(", method_name(info.config), ": ",
          info.n_interior, " interior / ", info.n_membrane, " membrane / ",
          info.n_outside, " outside of ", info.n_emitters, " emitters, ",
          round(info.runtime_s; digits = 2), "s)")
end
