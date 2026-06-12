"""
Internal: compute per-loop diagnostics for `loop_diagnostics.csv`.

Density at vertices is evaluated against the **originals-only** KDTree —
reflected/augmented points are used only for polygon construction, never
for diagnostics.
"""

const _DIAG_DENSITY_K = 128  # K used for vertex-density diagnostic
const _DIAG_DENSE_THRESH_FRAC = 0.8  # provisional heuristic_type threshold

function _heuristic_type(loop_id::Int, frac_in_fov::Float64,
                         frac_dense::Float64)::String
    loop_id == 1 && return "outer"
    frac_in_fov == 0.0 && return "reflection_noise"
    frac_in_fov < 1.0 && return "fov_crossing"
    frac_dense >= _DIAG_DENSE_THRESH_FRAC ? "interior_dense" : "interior_sparse"
end

function _compute_loop_diagnostics(
    loops::Vector{Vector{NTuple{2,Float64}}},
    x_um::Vector{Float64}, y_um::Vector{Float64},
    Xorig::Matrix{Float64},
    fov_um::NTuple{4,Float64},
    rho_thresh::Float64,
)
    fxmin, fxmax, fymin, fymax = fov_um
    n_orig = length(x_um)
    tree = NearestNeighbors.KDTree(Xorig)
    K = _DIAG_DENSITY_K
    inv_pi = 1 / π

    diags = Vector{LoopDiagnostic}(undef, length(loops))
    for (lid, verts) in enumerate(loops)
        nv = length(verts)
        area = abs(_polygon_area(verts))

        # n_emitters_inside (originals strictly inside)
        n_inside = 0
        @inbounds for i in 1:n_orig
            _point_in_polygon(x_um[i], y_um[i], verts) && (n_inside += 1)
        end

        # vertex density vs originals KDTree
        rhos = Vector{Float64}(undef, nv)
        in_fov_count = 0
        @inbounds for vi in 1:nv
            vx, vy = verts[vi]
            (fxmin <= vx <= fxmax && fymin <= vy <= fymax) && (in_fov_count += 1)
            _, dists = NearestNeighbors.knn(tree, [vx, vy], K, true)
            d = dists[end]
            rhos[vi] = (K - 1) * inv_pi / (d * d)
        end
        frac_dense = count(>=(rho_thresh), rhos) / nv
        med_rho = Statistics.median(rhos)
        frac_in_fov = in_fov_count / nv

        # v1 outer-polygon classifier: only loop_id == 1 participates in
        # inside_outer. For METHOD == "mask_carve", loop_id == 1 is still
        # marked `used_in_outer = true` because the alpha outer is the
        # source envelope used to build the carve — but the EFFECTIVE
        # classification boundary is the carve polygon, recorded in
        # `effective_outer.tsv` (not loop_id == 1 of polygon_loops.tsv).
        # Future "concave_refined" method may promote additional loops;
        # that decision will be made at result-build time.
        used_in_outer = (lid == 1)
        diags[lid] = LoopDiagnostic(
            lid, nv, area, n_inside, frac_in_fov, frac_dense, med_rho,
            used_in_outer,
            _heuristic_type(lid, frac_in_fov, frac_dense),
        )
    end
    return diags
end
