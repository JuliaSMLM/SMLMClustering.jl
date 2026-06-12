"""
Result records for `classify_emitters`.

`EdgeClassifyInfo <: SMLMData.AbstractSMLMInfo` is the sibling of `ClusterInfo` /
`ClusterStatisticsInfo`. Because emitters have no categorical class field (unlike
cluster `id`s), the canonical per-emitter answer lives here; `classify_emitters`
on an SMLD also mirrors the primary `class` vector into
`smld.metadata["edge_classify_class"]` for pipeline chaining.
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
    MaskCarveDiagnostic

Per-call diagnostic for `MaskCarveConfig`. `applied = false` means the carve fell
back to the v1 outer polygon (degenerate mask / empty intersection / polygonization
failure); the v1 polygon is used in that case.
"""
struct MaskCarveDiagnostic
    applied::Bool
    fallback_reason::String
    sigma_um::Float64
    k_noise::Float64
    pixel_um::Float64
    min_component_frac::Float64
    fill_hole_max_um2::Float64
    v1_polygon_area_um2::Float64
    carve_polygon_area_um2::Float64
    area_delta_um2::Float64
    v1_only_area_um2::Float64
    carve_only_area_um2::Float64
    med_v1_carve_distance_um::Float64
    p95_v1_carve_distance_um::Float64
    n_holes_filled::Int
    n_holes_preserved::Int
    n_carve_polygon_pts::Int
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

`inside_outer` is strictly **geometric** (containment in the alpha outer loop) with
`dist_to_outer_um` its distance (`NaN` when not inside). For `KdeValleyConfig` the
enclosure stage folds enclosed background into `class == :interior` while leaving
`inside_outer` geometric, so the enclosure-recovered set is exactly
`class == :interior && inside_outer == false` (with `dist_to_outer_um == NaN`).
Topological cell membership is `in_cell(info) == (class .!= :outside)`. Downstream
interior filters should read `class`, never `inside_outer`.

`config` holds the concrete config that ran (honest provenance);
`mask_carve_diagnostic` is populated only for `MaskCarveConfig`.
"""
struct EdgeClassifyInfo{C<:AbstractEdgeClassifyConfig} <: SMLMData.AbstractSMLMInfo
    n_emitters::Int
    class::Vector{Symbol}
    inside_outer::BitVector
    dist_to_outer_um::Vector{Float64}
    outer_polygon::Vector{NTuple{2,Float64}}
    loops::Vector{Vector{NTuple{2,Float64}}}
    loop_diagnostics::Vector{LoopDiagnostic}
    mask_carve_diagnostic::Union{Nothing,MaskCarveDiagnostic}
    config::C
    fov_um::NTuple{4,Float64}
    truncated_sides::NamedTuple{(:L, :R, :B, :T), NTuple{4,Bool}}
    n_reflected::Int
    runtime_s::Float64
    n_outside::Int
    n_membrane::Int
    n_interior::Int
end

"""
    in_cell(info) -> BitVector

Topological cell membership, `== (info.class .!= :outside)`. For
`KdeValleyConfig` this is a superset of `inside_outer` (includes enclosure-recovered
interior); for the polygon family it equals `inside_outer`.
"""
in_cell(info::EdgeClassifyInfo) = info.class .!= :outside

"""
    interior_fraction(info) -> Float64
"""
interior_fraction(info::EdgeClassifyInfo) = info.n_interior / max(info.n_emitters, 1)
