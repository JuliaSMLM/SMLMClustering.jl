"""
    compute_concavity_metric(result, x_um, y_um;
                             buffer_um = result.params_used.CONCAVITY_METRIC_BUFFER_NM/1000,
                             asym_R_nm = 1000.0,
                             asym_gate = 0.20,
                             rho_lo = 200.0,
                             intracellular_dense_threshold = 0.8)
        -> ConcavityMetricReport

Boundary-proximal concavity-error metric for the v1 outer-polygon
classifier. Identifies emitters that v1 calls `interior` but live in deep
concave bays the alpha-shape bridged across.

A v1 `interior` emitter is a **suspect** iff:

1. It lies within `buffer_um` of the outer polygon boundary (boundary-
   proximal — bays are at the membrane, not deep in the cell).
2. It is NOT inside any "intracellular void" loop, defined raw-metric-wise
   as a non-outer loop with `frac_in_fov == 1` AND
   `frac_dense >= intracellular_dense_threshold`. These are nuclei /
   sparse intracellular regions; they should not count as membrane
   concavity errors.
3. Its directional asymmetry at radius `asym_R_nm` is `>= asym_gate`
   (high asym → it lives near a real edge, contradicting v1's `interior`
   call).
4. Its local density `ρ_K(K=128)` (originals KDTree) is `<= rho_lo`
   (sparse local → not deep tissue).

Stratification: each suspect's nearest outer-polygon segment is
classified `interior_fov` if both endpoints are inside the FOV,
`fov_edge` otherwise. Approach B (vertex snap to asym ridge) is expected
to fix `interior_fov` suspects; Approach C (FOV-crossing loops as
auxiliary boundary) is expected to fix `fov_edge` suspects.

`buffer_um` defaults to `params_used.CONCAVITY_METRIC_BUFFER_NM/1000`.
"""
function compute_concavity_metric(
    result::EdgeClassificationResult,
    x_um::AbstractVector{<:Real},
    y_um::AbstractVector{<:Real};
    buffer_um::Float64 = result.params_used.CONCAVITY_METRIC_BUFFER_NM / 1000,
    asym_R_nm::Float64 = 1000.0,
    asym_gate::Float64 = 0.20,
    rho_lo::Float64 = 200.0,
    intracellular_dense_threshold::Float64 = 0.8,
)
    n = length(x_um)
    n == result.n_emitters ||
        throw(ArgumentError("x_um/y_um length must match result.n_emitters"))
    n == length(y_um) ||
        throw(ArgumentError("x_um and y_um must have equal length"))

    fxmin, fxmax, fymin, fymax = result.fov_um

    # Identify intracellular-void loops (raw-metric, not heuristic_type).
    # Note: result.loop_diagnostics[1] is the outer; never excluded.
    void_loop_ids = Int[]
    for d in result.loop_diagnostics
        d.loop_id == 1 && continue
        if d.frac_in_fov == 1.0 && d.frac_dense >= intracellular_dense_threshold
            push!(void_loop_ids, d.loop_id)
        end
    end
    void_polys = [result.loops[lid] for lid in void_loop_ids]

    # Originals KDTree for asym + ρ_K.
    Xorig = Matrix{Float64}(undef, 2, n)
    @inbounds for i in 1:n
        Xorig[1, i] = x_um[i]; Xorig[2, i] = y_um[i]
    end
    tree = NearestNeighbors.KDTree(Xorig)

    asym_R_um = asym_R_nm / 1000
    inv_pi = 1 / π

    n_interior = sum(==("interior"), result.class)
    n_eligible = 0
    suspects = Int[]
    suspect_is_fov_edge = Bool[]

    outer = result.outer_polygon
    no_outer = length(outer)

    @inbounds for i in 1:n
        result.class[i] == "interior" || continue
        d_outer = result.dist_to_outer_um[i]
        isnan(d_outer) && continue
        d_outer <= buffer_um || continue
        # Exclude intracellular voids
        in_void = false
        for vp in void_polys
            if _point_in_polygon(x_um[i], y_um[i], vp)
                in_void = true; break
            end
        end
        in_void && continue
        n_eligible += 1

        # Asymmetry at asym_R_nm.
        idxs = NearestNeighbors.inrange(tree, [x_um[i], y_um[i]], asym_R_um)
        m = length(idxs)
        m >= 4 || continue   # need a few neighbors for a meaningful centroid
        cx = 0.0; cy = 0.0
        for j in idxs
            cx += x_um[j]; cy += y_um[j]
        end
        cx /= m; cy /= m
        asym = hypot(cx - x_um[i], cy - y_um[i]) / asym_R_um
        asym >= asym_gate || continue

        # Local density ρ_K(K=128).
        _, dists = NearestNeighbors.knn(tree, [x_um[i], y_um[i]], 129, true)
        d128 = dists[end]
        rho = (128 - 1) * inv_pi / (d128 * d128)
        rho <= rho_lo || continue

        # Stratify by nearest outer segment.
        # Find nearest segment: argmin over i of dist to segment (i, i+1).
        best_seg = 1; best_d = Inf
        for s in 1:no_outer
            sj = s == no_outer ? 1 : s + 1
            ax_, ay_ = outer[s]; bx_, by_ = outer[sj]
            dx = bx_ - ax_; dy = by_ - ay_
            len2 = dx*dx + dy*dy
            if len2 == 0
                d = hypot(x_um[i] - ax_, y_um[i] - ay_)
            else
                t = ((x_um[i] - ax_) * dx + (y_um[i] - ay_) * dy) / len2
                t = clamp(t, 0.0, 1.0)
                cxs = ax_ + t * dx; cys = ay_ + t * dy
                d = hypot(x_um[i] - cxs, y_um[i] - cys)
            end
            if d < best_d; best_d = d; best_seg = s; end
        end
        sj = best_seg == no_outer ? 1 : best_seg + 1
        avx, avy = outer[best_seg]; bvx, bvy = outer[sj]
        a_in_fov = (fxmin <= avx <= fxmax) && (fymin <= avy <= fymax)
        b_in_fov = (fxmin <= bvx <= fxmax) && (fymin <= bvy <= fymax)
        is_fov_edge = !(a_in_fov && b_in_fov)

        push!(suspects, i)
        push!(suspect_is_fov_edge, is_fov_edge)
    end

    n_suspect = length(suspects)
    sx = Float64[x_um[i] for i in suspects]
    sy = Float64[y_um[i] for i in suspects]
    bv = BitVector(suspect_is_fov_edge)
    n_fov_edge = count(bv)
    n_int_fov = n_suspect - n_fov_edge

    return ConcavityMetricReport(
        buffer_um, asym_R_nm, asym_gate, rho_lo,
        n_interior, n_eligible, n_suspect, n_int_fov, n_fov_edge,
        sx, sy, bv,
    )
end
