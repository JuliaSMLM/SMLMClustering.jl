"""
Internal geometry helpers for `EdgeClassify`: FOV reflection, multi-K
density, alpha-shape boundary loops, polygon point-in / distance.
"""

# ---- FOV reflection ----------------------------------------------------------

function _truncated_sides(x_um, y_um, fov_um::NTuple{4,Float64}, tol_um::Float64)
    xmn, xmx = extrema(x_um); ymn, ymx = extrema(y_um)
    fxmin, fxmax, fymin, fymax = fov_um
    return (L = (xmn - fxmin) <= tol_um,
            R = (fxmax - xmx) <= tol_um,
            B = (ymn - fymin) <= tol_um,
            T = (fymax - ymx) <= tol_um)
end

function _reflect_emitters(x_um, y_um, fov_um::NTuple{4,Float64},
                           sides::NamedTuple, max_dist_um::Float64)
    fxmin, fxmax, fymin, fymax = fov_um
    n = length(x_um)
    refl_x = Float64[]; refl_y = Float64[]
    @inbounds for i in 1:n
        x, y = x_um[i], y_um[i]
        in_l = sides.L && (x - fxmin) < max_dist_um
        in_r = sides.R && (fxmax - x) < max_dist_um
        in_b = sides.B && (y - fymin) < max_dist_um
        in_t = sides.T && (fymax - y) < max_dist_um
        if in_l; push!(refl_x, 2*fxmin - x); push!(refl_y, y); end
        if in_r; push!(refl_x, 2*fxmax - x); push!(refl_y, y); end
        if in_b; push!(refl_x, x); push!(refl_y, 2*fymin - y); end
        if in_t; push!(refl_x, x); push!(refl_y, 2*fymax - y); end
        if in_l && in_b; push!(refl_x, 2*fxmin - x); push!(refl_y, 2*fymin - y); end
        if in_l && in_t; push!(refl_x, 2*fxmin - x); push!(refl_y, 2*fymax - y); end
        if in_r && in_b; push!(refl_x, 2*fxmax - x); push!(refl_y, 2*fymin - y); end
        if in_r && in_t; push!(refl_x, 2*fxmax - x); push!(refl_y, 2*fymax - y); end
    end
    xfull = vcat(x_um, refl_x); yfull = vcat(y_um, refl_y)
    return xfull, yfull, length(refl_x)
end

# ---- Multi-K density ---------------------------------------------------------

function _knn_K_density(X::AbstractMatrix{Float64}, K::Int, tree)
    n = size(X, 2)
    _, dists = NearestNeighbors.knn(tree, X, K + 1, true)
    rho = Vector{Float64}(undef, n)
    inv_pi = 1 / π
    @inbounds for i in 1:n
        d = dists[i][end]
        rho[i] = (K - 1) * inv_pi / (d * d)
    end
    return rho
end

# `k_list` is any iterable of integers (Tuple or Vector). An empty `k_list`
# disables the gate (returns all-true, no k-NN) — used by the kde_valley path.
# Each K is clamped to n_total-1 so small clouds don't throw on knn(K+1).
function _tissue_mask(X::Matrix{Float64}, k_list, rho_thresh::Float64)
    n_total = size(X, 2)
    (isempty(k_list) || n_total <= 1) && return trues(n_total)
    tree = NearestNeighbors.KDTree(X)
    mask = trues(n_total)
    for K in k_list
        Keff = min(Int(K), n_total - 1)
        Keff >= 1 || continue
        rho = _knn_K_density(X, Keff, tree)
        @inbounds for i in 1:n_total
            rho[i] >= rho_thresh || (mask[i] = false)
        end
    end
    return mask
end

# ---- Alpha-shape, boundary tracing, polygon geometry -------------------------

function _circumradius(p1, p2, p3)
    ax, ay = p1; bx, by = p2; cx, cy = p3
    d = 2 * (ax * (by - cy) + bx * (cy - ay) + cx * (ay - by))
    abs(d) < 1e-30 && return Inf
    a2 = ax^2 + ay^2; b2 = bx^2 + by^2; c2 = cx^2 + cy^2
    ux = (a2 * (by - cy) + b2 * (cy - ay) + c2 * (ay - by)) / d
    uy = (a2 * (cx - bx) + b2 * (ax - cx) + c2 * (bx - ax)) / d
    return hypot(ux - ax, uy - ay)
end

function _trace_boundary_loops(boundary_edges::Vector{Tuple{Int,Int}})
    adj = Dict{Int, Vector{Int}}()
    for (a, b) in boundary_edges
        push!(get!(adj, a, Int[]), b)
        push!(get!(adj, b, Int[]), a)
    end
    used = Set{Tuple{Int,Int}}()
    edge_key(a, b) = a < b ? (a, b) : (b, a)
    loops = Vector{Vector{Int}}()
    for (a0, b0) in boundary_edges
        edge_key(a0, b0) in used && continue
        loop = Int[a0, b0]
        push!(used, edge_key(a0, b0))
        prev, cur = a0, b0
        closed = false
        while true
            nbrs = get(adj, cur, Int[])
            nxt = 0
            for v in nbrs
                v == prev && continue
                edge_key(cur, v) in used && continue
                nxt = v; break
            end
            if nxt == 0 || nxt == a0
                if nxt == a0
                    push!(used, edge_key(cur, a0))
                    closed = true
                end
                break
            end
            push!(used, edge_key(cur, nxt))
            push!(loop, nxt)
            prev, cur = cur, nxt
        end
        # Keep only loops that actually closed back to the start; a closed ring
        # needs ≥3 distinct vertices (stored open, so a triangle has length 3).
        closed && length(loop) >= 3 && push!(loops, loop)
    end
    return loops
end

function _polygon_area(verts::Vector{NTuple{2,Float64}})
    n = length(verts); s = 0.0
    @inbounds for i in 1:n
        j = i == n ? 1 : i + 1
        s += verts[i][1] * verts[j][2] - verts[j][1] * verts[i][2]
    end
    return s / 2
end

# A polygon is usable (drawable / a valid Hopkins window) only with ≥3 DISTINCT
# vertices and a non-negligible area — Sutherland–Hodgman can emit duplicate
# on-edge vertices, so a raw vertex count alone is not a sufficient guard.
_valid_polygon(p::Vector{NTuple{2,Float64}}) =
    length(unique(p)) >= 3 && abs(_polygon_area(p)) > 1e-9

"""
    _alpha_shape_loops(X, alpha_um) -> Vector{Vector{NTuple{2,Float64}}}

Compute alpha-shape boundary loops on a 2×N point matrix `X`. Returns a
vector of closed-loop vertex sequences, sorted by `abs(polygon_area)`
descending. The first entry is the outer mosaic boundary; subsequent
entries are interior / hole / reflection-space loops.
"""
function _alpha_shape_loops(X::Matrix{Float64}, alpha_um::Float64)
    n = size(X, 2)
    points = [(X[1, i], X[2, i]) for i in 1:n]
    triangles = _delaunay_triangles(X)
    # Degenerate guard: <3 distinct points (or all-collinear) ⇒ no triangles ⇒
    # empty loops; the caller raises a clean error.
    isempty(triangles) && return Vector{Vector{NTuple{2,Float64}}}()
    edge_count = Dict{Tuple{Int,Int}, Int}()
    for (i, j, k) in triangles
        p1 = points[i]; p2 = points[j]; p3 = points[k]
        _circumradius(p1, p2, p3) <= alpha_um || continue
        for (a, b) in ((i, j), (j, k), (k, i))
            e = a < b ? (a, b) : (b, a)
            edge_count[e] = get(edge_count, e, 0) + 1
        end
    end
    boundary = Tuple{Int,Int}[(e[1], e[2]) for (e, c) in edge_count if c == 1]
    loop_idx = _trace_boundary_loops(boundary)
    polys = [[points[v] for v in idx] for idx in loop_idx]
    # Descending |area|, tie-broken by the lexicographically smallest vertex so
    # equal-area loops cannot reorder across runs / Julia versions (loops[1] must
    # be deterministic — the caller classifies against it).
    sort!(polys; by = p -> (-abs(_polygon_area(p)), minimum(p)), alg = MergeSort)
    return polys
end

# ---- FOV clipping (Sutherland–Hodgman) ---------------------------------------

# Intersection of segment a→b with the vertical line x = xc (called only when a, b
# straddle xc, so b[1] != a[1]); and with the horizontal line y = yc (b[2] != a[2]).
_seg_x(a::NTuple{2,Float64}, b::NTuple{2,Float64}, xc::Float64) =
    (xc, a[2] + (xc - a[1]) / (b[1] - a[1]) * (b[2] - a[2]))
_seg_y(a::NTuple{2,Float64}, b::NTuple{2,Float64}, yc::Float64) =
    (a[1] + (yc - a[2]) / (b[2] - a[2]) * (b[1] - a[1]), yc)

# One Sutherland–Hodgman pass: keep the `inside` half-plane, inserting the
# boundary intersection (`isect`) wherever the polygon edge crosses it.
function _clip_halfplane(input::Vector{NTuple{2,Float64}}, inside, isect)
    m = length(input)
    out = NTuple{2,Float64}[]
    m == 0 && return out
    @inbounds for i in 1:m
        cur = input[i]
        prev = input[i == 1 ? m : i - 1]
        cur_in = inside(cur); prev_in = inside(prev)
        if cur_in
            prev_in || push!(out, isect(prev, cur))
            push!(out, cur)
        elseif prev_in
            push!(out, isect(prev, cur))
        end
    end
    return out
end

"""
    _fov_clip(poly, fov_um) -> Vector{NTuple{2,Float64}}

Clip a closed polygon to the axis-aligned FOV rectangle
`(xmin, xmax, ymin, ymax)` via Sutherland–Hodgman. The FOV is convex, so the
result is the exact intersection `poly ∩ FOV`; vertices on an edge are kept.
Returns an empty vector when the intersection is empty.
"""
function _fov_clip(poly::Vector{NTuple{2,Float64}}, fov_um::NTuple{4,Float64})
    isempty(poly) && return poly
    fxmin, fxmax, fymin, fymax = fov_um
    out = poly
    out = _clip_halfplane(out, p -> p[1] >= fxmin, (a, b) -> _seg_x(a, b, fxmin))
    out = _clip_halfplane(out, p -> p[1] <= fxmax, (a, b) -> _seg_x(a, b, fxmax))
    out = _clip_halfplane(out, p -> p[2] >= fymin, (a, b) -> _seg_y(a, b, fymin))
    out = _clip_halfplane(out, p -> p[2] <= fymax, (a, b) -> _seg_y(a, b, fymax))
    return out
end

# ---- Perpendicular distance to polygon ---------------------------------------
# (`_point_in_polygon` is a shared helper in src/utils.jl, imported by this module.)

function _dist_to_polygon(qx::Float64, qy::Float64,
                          verts::Vector{NTuple{2,Float64}})
    n = length(verts); best = Inf
    @inbounds for i in 1:n
        j = i == n ? 1 : i + 1
        ax_, ay_ = verts[i]; bx_, by_ = verts[j]
        dx = bx_ - ax_; dy = by_ - ay_
        len2 = dx*dx + dy*dy
        if len2 == 0
            d = hypot(qx - ax_, qy - ay_)
        else
            t = ((qx - ax_) * dx + (qy - ay_) * dy) / len2
            t = clamp(t, 0.0, 1.0)
            cx = ax_ + t * dx; cy = ay_ + t * dy
            d = hypot(qx - cx, qy - cy)
        end
        d < best && (best = d)
    end
    return best
end

# ---- multi-cell pipeline helpers ---------------------------------------------

# Per-point k-th nearest-neighbor distance over a 2×N point set — the empirical
# local inter-point spacing (∝ 1/√ρ). `k` is clamped to n-1 for tiny clouds.
function _knn_distances(X::Matrix{Float64}, k::Int)
    n = size(X, 2)
    keff = min(k, n - 1)
    keff < 1 && return zeros(n)
    tree = NearestNeighbors.KDTree(X)
    _, dists = NearestNeighbors.knn(tree, X, keff + 1, true)
    return Float64[d[end] for d in dists]
end

# Middle of three values.
_median3(a, b, c) = a > b ? (b > c ? b : (a > c ? c : a)) : (a > c ? a : (b > c ? c : b))

# Multi-scale (per-triangle) alpha-shape loops. Each triangle's circumradius cap is
# `min(scale × median(k-NN of its 3 vertices), cap_um)` — the intersection of a LOCAL
# adaptive α (the scale×median term: SHRINKS where dense to carve real concavities,
# GROWS where sparse to bridge inter-clump gaps) with a CONSERVATIVE per-cell envelope
# `cap_um`. Because it is a min(), it can only remove triangles from the pure-local
# shape, never add bridges: the envelope rejects far-reaching low-density-noise
# protrusions (local α blows up at isolated noise but is capped), while a loose enough
# `cap_um` still lets a diffuse background coalesce. `kdist` is the per-point k-NN
# distance over the same point set `X`.
function _local_alpha_shape_loops(X::Matrix{Float64}, kdist::Vector{Float64},
                                  scale::Float64, cap_um::Float64)
    n = size(X, 2)
    points = [(X[1, i], X[2, i]) for i in 1:n]
    triangles = _delaunay_triangles(X)
    isempty(triangles) && return Vector{Vector{NTuple{2,Float64}}}()
    edge_count = Dict{Tuple{Int,Int}, Int}()
    for (i, j, k) in triangles
        a_local = min(scale * _median3(kdist[i], kdist[j], kdist[k]), cap_um)
        _circumradius(points[i], points[j], points[k]) <= a_local || continue
        for (a, b) in ((i, j), (j, k), (k, i))
            e = a < b ? (a, b) : (b, a)
            edge_count[e] = get(edge_count, e, 0) + 1
        end
    end
    boundary = Tuple{Int,Int}[(e[1], e[2]) for (e, c) in edge_count if c == 1]
    loop_idx = _trace_boundary_loops(boundary)
    polys = [[points[v] for v in idx] for idx in loop_idx]
    sort!(polys; by = p -> (-abs(_polygon_area(p)), minimum(p)), alg = MergeSort)
    return polys
end

# Relative-density gate: of tissue indices `idx`, keep those whose neighbor count
# within `r_um` is at least `frac × median(count)`. Removes isolated outlier
# whiskers without erasing genuine low-density cells (the cutoff is relative to the
# cell's own median, so it self-scales with density). `frac <= 0` disables it.
function _relative_core_filter(x::Vector{Float64}, y::Vector{Float64},
                               idx::Vector{Int}, r_um::Float64, frac::Float64)
    (frac <= 0 || length(idx) < 3) && return idx
    m = length(idx)
    X = Matrix{Float64}(undef, 2, m)
    @inbounds for k in 1:m
        X[1, k] = x[idx[k]]; X[2, k] = y[idx[k]]
    end
    tree = NearestNeighbors.KDTree(X)
    cnt = Vector{Int}(undef, m)
    @inbounds for k in 1:m
        cnt[k] = length(NearestNeighbors.inrange(tree, view(X, :, k), r_um)) - 1
    end
    thr = frac * median(cnt)
    return Int[idx[k] for k in 1:m if cnt[k] >= thr]
end

# A boundary segment a→b is a field-of-view CUT (not a real membrane) when BOTH
# endpoints sit within `tol` of the same FOV side AND that side is actually
# truncated (`sides`). Gating on `sides` prevents a real cell edge that merely runs
# near a non-truncated FOV border from being mistaken for a cut. `<= tol` so an
# exact-edge tolerance of 0 still matches points lying on the edge.
function _seg_on_fov_edge(a::NTuple{2,Float64}, b::NTuple{2,Float64},
                          fov::NTuple{4,Float64}, tol::Float64,
                          sides::NamedTuple)
    fxmin, fxmax, fymin, fymax = fov
    (sides.L && abs(a[1] - fxmin) <= tol && abs(b[1] - fxmin) <= tol) && return true
    (sides.R && abs(a[1] - fxmax) <= tol && abs(b[1] - fxmax) <= tol) && return true
    (sides.B && abs(a[2] - fymin) <= tol && abs(b[2] - fymin) <= tol) && return true
    (sides.T && abs(a[2] - fymax) <= tol && abs(b[2] - fymax) <= tol) && return true
    return false
end

# Distance from a query point to segment a→b.
function _dist_point_seg(qx::Float64, qy::Float64,
                         a::NTuple{2,Float64}, b::NTuple{2,Float64})
    dx = b[1] - a[1]; dy = b[2] - a[2]
    len2 = dx * dx + dy * dy
    len2 == 0 && return hypot(qx - a[1], qy - a[2])
    t = clamp(((qx - a[1]) * dx + (qy - a[2]) * dy) / len2, 0.0, 1.0)
    return hypot(qx - (a[1] + t * dx), qy - (a[2] + t * dy))
end

# Distance from (qx,qy) to the nearest segment of a single ring that is NOT on a
# truncated FOV edge. Returns Inf when every segment lies on a FOV edge (no real
# membrane there).
function _dist_to_ring_excl_fov(qx::Float64, qy::Float64,
                                ring::Vector{NTuple{2,Float64}},
                                fov::NTuple{4,Float64}, tol::Float64,
                                sides::NamedTuple)
    m = length(ring)
    m < 2 && return Inf
    best = Inf
    @inbounds for i in 1:m
        a = ring[i]; b = ring[i == m ? 1 : i + 1]
        _seg_on_fov_edge(a, b, fov, tol, sides) && continue
        d = _dist_point_seg(qx, qy, a, b)
        d < best && (best = d)
    end
    return best
end
