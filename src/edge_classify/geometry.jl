"""
Internal geometry helpers for `EdgeClassify`: FOV reflection, multi-K
density, alpha-shape boundary loops, polygon point-in / distance.
"""

# ---- FOV reflection ----------------------------------------------------------

function _truncated_sides(x_um, y_um, fov_um::NTuple{4,Float64}, tol_um::Float64)
    xmn, xmx = extrema(x_um); ymn, ymx = extrema(y_um)
    fxmin, fxmax, fymin, fymax = fov_um
    return (L = (xmn - fxmin) < tol_um,
            R = (fxmax - xmx) < tol_um,
            B = (ymn - fymin) < tol_um,
            T = (fymax - ymx) < tol_um)
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
        while true
            nbrs = get(adj, cur, Int[])
            nxt = 0
            for v in nbrs
                v == prev && continue
                edge_key(cur, v) in used && continue
                nxt = v; break
            end
            if nxt == 0 || nxt == a0
                nxt == a0 && push!(used, edge_key(cur, a0))
                break
            end
            push!(used, edge_key(cur, nxt))
            push!(loop, nxt)
            prev, cur = cur, nxt
        end
        length(loop) >= 4 && push!(loops, loop)
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
    tri = DelaunayTriangulation.triangulate(points)
    edge_count = Dict{Tuple{Int,Int}, Int}()
    for T in DelaunayTriangulation.each_solid_triangle(tri)
        i, j, k = DelaunayTriangulation.triangle_vertices(T)
        p1 = DelaunayTriangulation.get_point(tri, i)
        p2 = DelaunayTriangulation.get_point(tri, j)
        p3 = DelaunayTriangulation.get_point(tri, k)
        _circumradius(p1, p2, p3) <= alpha_um || continue
        for (a, b) in ((i, j), (j, k), (k, i))
            e = a < b ? (a, b) : (b, a)
            edge_count[e] = get(edge_count, e, 0) + 1
        end
    end
    boundary = Tuple{Int,Int}[(e[1], e[2]) for (e, c) in edge_count if c == 1]
    loop_idx = _trace_boundary_loops(boundary)
    polys = [[points[v] for v in idx] for idx in loop_idx]
    sort!(polys; by = p -> -abs(_polygon_area(p)), alg = MergeSort)  # stable → reproducible loop order
    return polys
end

# ---- Point-in-polygon, perpendicular distance --------------------------------

function _point_in_polygon(qx::Float64, qy::Float64,
                           verts::Vector{NTuple{2,Float64}})
    n = length(verts); inside = false; j = n
    @inbounds for i in 1:n
        ix, iy = verts[i]; jx, jy = verts[j]
        if (iy > qy) != (jy > qy)
            xint = (jx - ix) * (qy - iy) / (jy - iy) + ix
            qx < xint && (inside = !inside)
        end
        j = i
    end
    return inside
end

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
