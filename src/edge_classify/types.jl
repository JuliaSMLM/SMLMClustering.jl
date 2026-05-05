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
end

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
    heuristic_type::String
end

"""
    EdgeClassificationResult

Result of `classify_emitters`. Class labels partition the input set:
`"outside" ∪ "membrane" ∪ "interior" == 1:n_emitters`.
"""
struct EdgeClassificationResult
    n_emitters::Int
    class::Vector{String}
    inside_outer::BitVector
    dist_to_outer_um::Vector{Float64}             # NaN where inside_outer == false
    outer_polygon::Vector{NTuple{2,Float64}}      # closed-loop vertices
    loops::Vector{Vector{NTuple{2,Float64}}}      # all loops, loop_id == index
    loop_diagnostics::Vector{LoopDiagnostic}
    params_used::EdgeClassifyParams
    fov_um::NTuple{4,Float64}
    truncated_sides::NamedTuple{(:L, :R, :B, :T), NTuple{4,Bool}}
    n_reflected::Int
    runtime_s::Float64
end
