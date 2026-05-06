"""
    EdgeClassifyParams

Parameters for the v1 outer-polygon edge/membrane/interior classifier.
All keys uppercase to match the documented `params.toml` / `params.json`
convention. Defaults are provisional and may move; callers pinning a
specific set should record `params_used` from the result.
"""
Base.@kwdef struct EdgeClassifyParams
    K_LIST::Vector{Int}            = [16, 128]
    RHO_K_THRESH::Float64          = 200.0
    ALPHA_NM::Float64              = 300.0
    REFLECT_RADIUS_NM::Float64     = 1500.0
    MEMBRANE_NM::Float64           = 100.0
    FOV_TRUNC_TOL_NM::Float64      = 150.0
    # Method selector. Default "outer_polygon" preserves v1 behavior bit-for-bit.
    # "grid_hybrid" preserves the v1 outside/interior topology and promotes
    # only v1-interior emitters that lie on a local density-grid boundary near
    # the v1 outer polygon.
    METHOD::String                 = "outer_polygon"
    GRID_PX_NM::Float64            = 50.0
    GRID_SMOOTH_NM::Float64        = 80.0
    GRID_MASK_Q::Float64           = 0.03
    GRID_MASK_PEAK_FRAC::Float64   = 0.26
    GRID_OUTER_BUFFER_NM::Float64  = 800.0
    # Concavity-evaluation buffer. Interior emitters within this distance of the
    # outer polygon are eligible to be flagged as "suspect" by the concavity
    # metric. Used by `compute_concavity_metric` only; does not affect class.
    CONCAVITY_METRIC_BUFFER_NM::Float64 = 2000.0
    # mask_carve method (provisional/opt-in). Carves v1 outer polygon inward
    # using a heterogeneity-robust density mask. Never expands v1 outward.
    # Defaults from dev/scripts/mask_contour_v3.jl synthetic-best (no real-cell
    # tuning). See docs/src/edge_classify_interface_v1.md for the carve-only
    # limitation.
    MASK_CARVE_SIGMA_UM::Float64           = 0.080
    MASK_CARVE_K_NOISE::Float64            = 3.0
    MASK_CARVE_PIXEL_UM::Float64           = 0.040
    MASK_CARVE_MIN_COMPONENT_FRAC::Float64 = 0.05
    MASK_CARVE_FILL_HOLE_MAX_UM2::Float64  = 0.5
end

const _METHOD_OUTER_POLYGON   = "outer_polygon"
const _METHOD_CONCAVE_REFINED = "concave_refined"
const _METHOD_GRID_HYBRID     = "grid_hybrid"
const _METHOD_MASK_CARVE      = "mask_carve"
const _VALID_METHODS = (_METHOD_OUTER_POLYGON, _METHOD_GRID_HYBRID,
                        _METHOD_CONCAVE_REFINED, _METHOD_MASK_CARVE)

"""
    LoopDiagnostic

Per-loop diagnostic record matching `loop_diagnostics.csv` column order.
`heuristic_type` column name is stable; values/thresholds are provisional.
"""
struct LoopDiagnostic
    loop_id::Int
    vertex_count::Int
    area_um2::Float64
    n_emitters_inside::Int
    frac_in_fov::Float64
    frac_dense::Float64
    median_rhoK::Float64
    used_in_outer::Bool   # schema_version 2: true iff this loop participates in the inside_outer decision
    heuristic_type::String
end

"""
    ConcavityMetricReport

Boundary-proximal concavity-error metric for the v1 outer-polygon
classifier. A "suspect" is a v1 `interior` emitter that lives near the
outer polygon, has high directional asymmetry at long radius, and low
local density — i.e., it sits in a deep concave bay the alpha-shape
bridged across.

Suspects are stratified by whether the nearest outer-polygon segment is
inside the FOV (interior-FOV concavity, fixable by Approach B) or
straddles the FOV edge (FOV-edge concavity, fixable by Approach C).

Emitters inside diagnostic interior-dense loops (real intracellular
voids) are excluded so nuclei are not counted as membrane concavity
errors.
"""
struct ConcavityMetricReport
    buffer_um::Float64
    asym_R_nm::Float64
    asym_gate::Float64
    rho_lo::Float64
    n_interior::Int
    n_eligible::Int                       # within buffer of outer, not inside intracellular void
    n_suspect::Int
    n_suspect_interior_fov::Int
    n_suspect_fov_edge::Int
    suspect_x_um::Vector{Float64}
    suspect_y_um::Vector{Float64}
    suspect_is_fov_edge::BitVector
end

"""
    MaskCarveDiagnostic

Per-call diagnostic for METHOD = "mask_carve". `applied = false` indicates
the carve was attempted but fell back to v1 outer polygon (degenerate D
mask, empty intersection, or polygonization failure); the v1 classifier
output is preserved in that case.

Areas measured by raster integration over the FOV at MASK_CARVE_PIXEL_UM
pixel pitch. `carve_only_area_um2` should be ≈ 0 by construction (modulo
rasterization roundoff).
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
    area_delta_um2::Float64                # carve − v1
    v1_only_area_um2::Float64              # v1 ∧ ¬carve
    carve_only_area_um2::Float64           # carve ∧ ¬v1 (≈0 invariant)
    med_v1_carve_distance_um::Float64
    p95_v1_carve_distance_um::Float64
    n_holes_filled::Int
    n_holes_preserved::Int
    n_carve_polygon_pts::Int
end

"""
    EdgeClassificationResult

Result of `classify_emitters`. Class labels partition the input set:
`"outside" ∪ "membrane" ∪ "interior" == 1:n_emitters`.

For `METHOD == "mask_carve"`, `outer_polygon` is the **effective**
classification polygon (the carve), while `loops[1]` retains the
alpha-shape outer loop as provenance. For all other methods,
`outer_polygon == loops[1]`.
"""
struct EdgeClassificationResult
    n_emitters::Int
    class::Vector{String}
    inside_outer::BitVector
    dist_to_outer_um::Vector{Float64}             # NaN where inside_outer == false
    outer_polygon::Vector{NTuple{2,Float64}}      # closed-loop vertices (effective)
    loops::Vector{Vector{NTuple{2,Float64}}}      # all loops, loop_id == index
    loop_diagnostics::Vector{LoopDiagnostic}
    params_used::EdgeClassifyParams
    fov_um::NTuple{4,Float64}
    truncated_sides::NamedTuple{(:L, :R, :B, :T), NTuple{4,Bool}}
    n_reflected::Int
    runtime_s::Float64
    mask_carve_diagnostic::Union{Nothing, MaskCarveDiagnostic}
end
