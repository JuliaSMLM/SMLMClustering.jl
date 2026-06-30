"""
    compute_concavity_metric(info::EdgeClassifyInfo, x_um, y_um;
                             buffer_um = 2.0, asym_R_nm = 1000.0, asym_gate = 0.20,
                             rho_lo = 200.0, intracellular_dense_threshold = 0.8)
        -> ConcavityMetricReport

Boundary-proximal concavity-error metric for the outer-polygon classifier.
Identifies emitters classified `:interior` that live in deep concave bays the
alpha-shape bridged across.

A `:interior` emitter is a **suspect** iff: (1) it lies within `buffer_um` of the
outer polygon boundary; (2) it is not inside an "intracellular void" loop
(non-outer loop with `frac_in_fov == 1` and `frac_dense >= intracellular_dense_threshold`);
(3) its directional asymmetry at `asym_R_nm` is `>= asym_gate`; (4) its local
density ρ_K(K=128) is `<= rho_lo`.

`buffer_um` is a diagnostic parameter (default 2.0 µm), not part of the classifier
config.
"""
function compute_concavity_metric(
    info::EdgeClassifyInfo,
    x_um::AbstractVector{<:Real},
    y_um::AbstractVector{<:Real};
    buffer_um::Float64 = 2.0,
    asym_R_nm::Float64 = 1000.0,
    asym_gate::Float64 = 0.20,
    rho_lo::Float64 = 200.0,
    intracellular_dense_threshold::Float64 = 0.8,
)
    n = length(x_um)
    n == info.n_emitters ||
        throw(ArgumentError("x_um/y_um length must match info.n_emitters"))
    n == length(y_um) ||
        throw(ArgumentError("x_um and y_um must have equal length"))

    fxmin, fxmax, fymin, fymax = info.fov_um

    # Intracellular-void loops (raw-metric). loop_diagnostics[1] is the outer.
    void_loop_ids = Int[]
    for d in info.loop_diagnostics
        d.loop_id == 1 && continue
        if d.frac_in_fov == 1.0 && d.frac_dense >= intracellular_dense_threshold
            push!(void_loop_ids, d.loop_id)
        end
    end
    void_polys = [info.loops[lid] for lid in void_loop_ids]

    Xorig = Matrix{Float64}(undef, 2, n)
    @inbounds for i in 1:n
        Xorig[1, i] = x_um[i]; Xorig[2, i] = y_um[i]
    end
    tree = NearestNeighbors.KDTree(Xorig)

    K_density = 128
    rho_all = _knn_K_density(Xorig, K_density, tree)

    asym_R_um = asym_R_nm / 1000

    n_interior = count(==(:interior), info.class)
    n_eligible = 0
    suspects = Int[]
    suspect_is_fov_edge = Bool[]

    # The classification boundary (reflected loop), NOT info.outer_polygon — the
    # latter is now the FOV-clipped published footprint, while dist_to_outer_um is
    # measured against this loop, so the two must stay paired here.
    outer = info.loops[1]
    no_outer = length(outer)

    @inbounds for i in 1:n
        info.class[i] == :interior || continue
        d_outer = info.dist_to_outer_um[i]
        isnan(d_outer) && continue
        d_outer <= buffer_um || continue
        in_void = false
        for vp in void_polys
            if _point_in_polygon(x_um[i], y_um[i], vp)
                in_void = true; break
            end
        end
        in_void && continue
        n_eligible += 1

        idxs = NearestNeighbors.inrange(tree, [x_um[i], y_um[i]], asym_R_um)
        m = length(idxs)
        m >= 4 || continue
        cx = 0.0; cy = 0.0
        for j in idxs
            cx += x_um[j]; cy += y_um[j]
        end
        cx /= m; cy /= m
        asym = hypot(cx - x_um[i], cy - y_um[i]) / asym_R_um
        asym >= asym_gate || continue

        rho_all[i] <= rho_lo || continue

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
