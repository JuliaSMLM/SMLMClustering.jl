"""
Pure-Julia 2D Delaunay triangulation (Bowyer–Watson incremental) for
`EdgeClassify`'s alpha-shape. Replaces the heavier `DelaunayTriangulation.jl`
engine in the alpha-shape path only.

Algorithm:
  1. Deduplicate the input points (exact `(x, y)` equality) so degenerate
     near/exact-coincident clouds cannot wedge the incremental insert.
  2. Build an enclosing super-triangle far outside the point bounding box.
  3. Order the insertion by a Hilbert space-filling curve so consecutive
     inserts are spatially local (keeps the point-location walk short).
  4. Insert each point: locate the containing triangle by *walking* from the
     last-inserted triangle (using `orient2`), flood-fill the cavity of all
     triangles whose circumcircle contains the point (`incircle`), delete the
     cavity, and re-triangulate the star-shaped hole by joining the new point
     to every cavity boundary edge.
  5. Drop every triangle incident to a super-triangle vertex.

Robustness comes entirely from `AdaptivePredicates`' exact `orient2` /
`incircle`; no float predicate is hand-rolled.
"""

# ---- triangle / adjacency storage -------------------------------------------
#
# Triangles are stored in flat arrays indexed by a triangle id `t`:
#   tri_v[t]   = (i, j, k)  vertex ids (into the working point list), CCW.
#   tri_adj[t] = (a, b, c)  neighbour triangle id across the edge OPPOSITE
#                           vertex 1, 2, 3 respectively (0 = no neighbour).
# Edge opposite vertex slot 1 is (j,k); slot 2 is (k,i); slot 3 is (i,j).
# `tri_dead[t]` marks deleted triangles (slots are recycled via `free`).

struct _DelaunayState
    px::Vector{Float64}        # working point x (incl. 3 super-triangle pts)
    py::Vector{Float64}        # working point y
    tri_v::Vector{NTuple{3,Int}}
    tri_adj::Vector{NTuple{3,Int}}
    tri_dead::Vector{Bool}
    free::Vector{Int}          # recycled triangle slots
end

@inline _pt(s::_DelaunayState, i::Int) = (s.px[i], s.py[i])

@inline function _orient(s::_DelaunayState, a::Int, b::Int, c::Int)
    return AdaptivePredicates.orient2(_pt(s, a), _pt(s, b), _pt(s, c))
end

@inline function _in_circle(s::_DelaunayState, a::Int, b::Int, c::Int, d::Int)
    return AdaptivePredicates.incircle(_pt(s, a), _pt(s, b), _pt(s, c), _pt(s, d))
end

function _new_triangle!(s::_DelaunayState, v::NTuple{3,Int}, adj::NTuple{3,Int})
    if isempty(s.free)
        push!(s.tri_v, v); push!(s.tri_adj, adj); push!(s.tri_dead, false)
        return length(s.tri_v)
    else
        t = pop!(s.free)
        s.tri_v[t] = v; s.tri_adj[t] = adj; s.tri_dead[t] = false
        return t
    end
end

@inline function _kill_triangle!(s::_DelaunayState, t::Int)
    s.tri_dead[t] = true
    push!(s.free, t)
end

# Set the neighbour of triangle `t` across the edge OPPOSITE its vertex slot
# `slot` (1,2,3) to `nbr`.
@inline function _set_adj!(s::_DelaunayState, t::Int, slot::Int, nbr::Int)
    a = s.tri_adj[t]
    s.tri_adj[t] = slot == 1 ? (nbr, a[2], a[3]) :
                   slot == 2 ? (a[1], nbr, a[3]) :
                               (a[1], a[2], nbr)
    return nothing
end

# ---- Hilbert spatial sort ----------------------------------------------------
#
# Map each point to a 16-bit-per-axis Hilbert index and sort by it. This gives
# strong insertion locality so the walk is O(1) amortised.

function _hilbert_d2xy_inv(order::Int, x0::UInt32, y0::UInt32)
    # Convert (x,y) -> Hilbert distance d for a 2^order grid.
    x = x0; y = y0
    rx = UInt32(0); ry = UInt32(0)
    d = UInt64(0)
    s = UInt32(1) << (order - 1)
    while s > 0
        rx = (x & s) > 0 ? UInt32(1) : UInt32(0)
        ry = (y & s) > 0 ? UInt32(1) : UInt32(0)
        d += UInt64(s) * UInt64(s) * UInt64((UInt32(3) * rx) ⊻ ry)
        # rotate
        if ry == 0
            if rx == 1
                x = s - 1 - x
                y = s - 1 - y
            end
            x, y = y, x
        end
        s >>= 1
    end
    return d
end

function _hilbert_order(px::AbstractVector{Float64}, py::AbstractVector{Float64},
                        ids::Vector{Int})
    n = length(ids)
    n == 0 && return ids
    xmn = Inf; xmx = -Inf; ymn = Inf; ymx = -Inf
    @inbounds for i in ids
        x = px[i]; y = py[i]
        x < xmn && (xmn = x); x > xmx && (xmx = x)
        y < ymn && (ymn = y); y > ymx && (ymx = y)
    end
    order = 16
    side = (UInt32(1) << order) - UInt32(1)
    sx = xmx > xmn ? Float64(side) / (xmx - xmn) : 0.0
    sy = ymx > ymn ? Float64(side) / (ymx - ymn) : 0.0
    keys = Vector{UInt64}(undef, n)
    @inbounds for (q, i) in enumerate(ids)
        # clamp guards a degenerate/huge span (sx→Inf) from NaN'ing the cast;
        # the Hilbert key only sets insertion ORDER, never the (unique) output.
        g = (px[i] - xmn) * sx
        gx = isfinite(g) ? UInt32(clamp(round(g), 0.0, Float64(side))) : UInt32(0)
        g = (py[i] - ymn) * sy
        gy = isfinite(g) ? UInt32(clamp(round(g), 0.0, Float64(side))) : UInt32(0)
        keys[q] = _hilbert_d2xy_inv(order, gx, gy)
    end
    perm = sortperm(keys)
    return ids[perm]
end

# ---- point location (walk) ---------------------------------------------------
#
# Walk from triangle `start` toward the triangle that contains point `p`
# (vertex id). At each triangle, if `p` is to the right of any directed edge
# (orient < 0 for a CCW triangle), step to the neighbour across that edge.
# Returns the containing triangle id. Falls back to a linear scan if the walk
# fails to terminate (defensive; should not happen with exact predicates).

function _locate(s::_DelaunayState, p::Int, start::Int)
    t = start
    (t == 0 || s.tri_dead[t]) && (t = _first_live(s))
    maxsteps = 2 * length(s.tri_v) + 16
    prev = 0
    for _ in 1:maxsteps
        v = s.tri_v[t]; a = s.tri_adj[t]
        # edge opposite slot1 = (v2,v3); slot2 = (v3,v1); slot3 = (v1,v2)
        o1 = _orient(s, v[2], v[3], p)
        if o1 < 0 && a[1] != prev && a[1] != 0
            prev = t; t = a[1]; continue
        end
        o2 = _orient(s, v[3], v[1], p)
        if o2 < 0 && a[2] != prev && a[2] != 0
            prev = t; t = a[2]; continue
        end
        o3 = _orient(s, v[1], v[2], p)
        if o3 < 0 && a[3] != prev && a[3] != 0
            prev = t; t = a[3]; continue
        end
        # No edge strictly excludes p (or its neighbour is the one we came
        # from / a boundary): we are in (or on) this triangle.
        if o1 >= 0 && o2 >= 0 && o3 >= 0
            return t
        end
        # p is right of an edge but that step is blocked (boundary). Try any
        # remaining unvisited neighbour that excludes p.
        if o1 < 0 && a[1] != 0 && a[1] != prev
            prev = t; t = a[1]; continue
        elseif o2 < 0 && a[2] != 0 && a[2] != prev
            prev = t; t = a[2]; continue
        elseif o3 < 0 && a[3] != 0 && a[3] != prev
            prev = t; t = a[3]; continue
        end
        # Stuck; break to linear fallback.
        break
    end
    # Linear fallback: find any triangle that contains p.
    @inbounds for tt in 1:length(s.tri_v)
        s.tri_dead[tt] && continue
        v = s.tri_v[tt]
        if _orient(s, v[2], v[3], p) >= 0 &&
           _orient(s, v[3], v[1], p) >= 0 &&
           _orient(s, v[1], v[2], p) >= 0
            return tt
        end
    end
    # The linear scan above must find a container (p lies inside the super-
    # triangle). Reaching here means a triangulation invariant broke — fail
    # loudly rather than return a non-containing triangle and corrupt the cavity.
    error("_delaunay: point location failed (triangulation invariant violated)")
end

@inline function _first_live(s::_DelaunayState)
    @inbounds for t in 1:length(s.tri_v)
        s.tri_dead[t] || return t
    end
    return 0
end

# ---- insertion (cavity flood-fill + re-triangulation) ------------------------
#
# Boundary edges of the cavity are collected as a star-shaped polygon around p.
# Each boundary edge stores (eu, ev, outnbr) where (eu -> ev) is CCW as seen
# from inside the cavity and `outnbr` is the LIVE triangle just outside that
# edge (0 if the edge is on the super-triangle hull). We then fan-triangulate:
# new triangle (p, eu, ev) for every boundary edge, re-stitching adjacency with
# `outnbr` (matched by edge vertices, not by stale ids) and between fan spokes.

function _insert_point!(s::_DelaunayState, p::Int, start::Int,
                        cavity::Vector{Int}, bnd_u::Vector{Int},
                        bnd_v::Vector{Int}, bnd_out::Vector{Int})
    t0 = _locate(s, p, start)
    t0 == 0 && return start

    # Flood-fill the cavity: all triangles whose circumcircle contains p.
    empty!(cavity)
    empty!(bnd_u); empty!(bnd_v); empty!(bnd_out)
    push!(cavity, t0)
    s.tri_dead[t0] = true            # mark as "in cavity" (temporarily)
    head = 1
    while head <= length(cavity)
        t = cavity[head]; head += 1
        v = s.tri_v[t]; a = s.tri_adj[t]
        # For each edge (opposite slot), test the neighbour.
        for slot in 1:3
            nbr = slot == 1 ? a[1] : slot == 2 ? a[2] : a[3]
            # the edge for this slot (CCW boundary order of the cavity)
            eu, ev = slot == 1 ? (v[2], v[3]) :
                     slot == 2 ? (v[3], v[1]) :
                                 (v[1], v[2])
            if nbr == 0
                # super-triangle outer boundary edge: always a cavity boundary
                push!(bnd_u, eu); push!(bnd_v, ev); push!(bnd_out, 0)
            elseif s.tri_dead[nbr]
                # already in cavity (marked dead during this flood). Shared
                # interior edge — not a cavity boundary; skip.
                continue
            else
                nv = s.tri_v[nbr]
                if _in_circle(s, nv[1], nv[2], nv[3], p) > 0
                    push!(cavity, nbr)
                    s.tri_dead[nbr] = true
                else
                    push!(bnd_u, eu); push!(bnd_v, ev); push!(bnd_out, nbr)
                end
            end
        end
    end

    nb = length(bnd_u)
    # Recycle the dead cavity slots. Their adjacency is no longer read: outside
    # neighbours are repointed by matching edge vertices (`_repoint_edge!`), not
    # by referencing the dead cavity ids.
    for t in cavity
        push!(s.free, t)
    end

    # Fan-triangulate: new triangle (p, u, w) per boundary edge (u -> w).
    newtris = Vector{Int}(undef, nb)
    @inbounds for e in 1:nb
        u = bnd_u[e]; w = bnd_v[e]; outnbr = bnd_out[e]
        # slots: 1=p, 2=u, 3=w. Edge opposite slot1 is (u,w) -> outnbr.
        t = _new_triangle!(s, (p, u, w), (outnbr, 0, 0))
        newtris[e] = t
        # Repoint outnbr's edge (w,u) from the (dead) cavity tri to `t`.
        _repoint_edge!(s, outnbr, w, u, t)
    end

    # Stitch adjacency between fan triangles sharing a spoke edge (p, vertex).
    _stitch_fan!(s, newtris, bnd_u, bnd_v)

    return newtris[1]
end

# In LIVE triangle `t`, find the edge equal to the undirected pair {a,b} and set
# that edge's neighbour to `newnbr`. The edge opposite slot `k` is the pair of
# the OTHER two vertices. Robust: matches by vertex identity, not by stale ids.
@inline function _repoint_edge!(s::_DelaunayState, t::Int, a::Int, b::Int, newnbr::Int)
    t == 0 && return nothing
    v = s.tri_v[t]
    # edge opposite slot1 = {v2,v3}; slot2 = {v3,v1}; slot3 = {v1,v2}
    if (v[2] == a && v[3] == b) || (v[2] == b && v[3] == a)
        _set_adj!(s, t, 1, newnbr)
    elseif (v[3] == a && v[1] == b) || (v[3] == b && v[1] == a)
        _set_adj!(s, t, 2, newnbr)
    elseif (v[1] == a && v[2] == b) || (v[1] == b && v[2] == a)
        _set_adj!(s, t, 3, newnbr)
    end
    return nothing
end

function _stitch_fan!(s::_DelaunayState, newtris::Vector{Int},
                      bnd_u::Vector{Int}, bnd_v::Vector{Int})
    nb = length(newtris)
    # Each fan triangle (p,u,w) owns spoke edges (p,u) [slot3] and (p,w) [slot2].
    # A spoke endpoint vertex appears in exactly two fan triangles; pair them.
    spokes = Dict{Int,Tuple{Int,Int}}()  # vertex -> (triangle, slot), first seen
    @inbounds for e in 1:nb
        t = newtris[e]; u = bnd_u[e]; w = bnd_v[e]
        # spoke (p,u) lives in slot3; spoke (p,w) lives in slot2
        for (vert, slot) in ((u, 3), (w, 2))
            if haskey(spokes, vert)
                t2, slot2 = spokes[vert]
                _set_adj!(s, t, slot, t2)
                _set_adj!(s, t2, slot2, t)
                delete!(spokes, vert)
            else
                spokes[vert] = (t, slot)
            end
        end
    end
    return nothing
end

# ---- driver ------------------------------------------------------------------

"""
    _delaunay_triangles(X::Matrix{Float64}) -> Vector{NTuple{3,Int}}

Bowyer–Watson Delaunay triangulation of the columns of `X` (2×N). Returns
triangles as 1-based vertex-index triples into the *columns of X* (CCW). Exact
duplicate points are collapsed (the first column index is used). Returns an
empty vector if fewer than 3 distinct points (or all points collinear).
"""
function _delaunay_triangles(X::Matrix{Float64})
    n = size(X, 2)
    n < 3 && return NTuple{3,Int}[]

    # Deduplicate exact (x,y) — map each distinct point to its first column id.
    coord_to_id = Dict{Tuple{Float64,Float64},Int}()
    uniq_ids = Int[]                     # original column id of each distinct pt
    sizehint!(coord_to_id, n)
    @inbounds for i in 1:n
        xi = X[1, i]; yi = X[2, i]
        (isfinite(xi) & isfinite(yi)) ||
            throw(ArgumentError("_delaunay_triangles: non-finite coordinate at column $i"))
        # Canonicalize -0.0 → 0.0: isequal(0.0,-0.0) is false, so without this a
        # signed-zero would survive dedup as a phantom coincident vertex and
        # corrupt the triangulation. (-0.0 == 0.0 is true, so iszero catches it.)
        key = (iszero(xi) ? 0.0 : xi, iszero(yi) ? 0.0 : yi)
        if !haskey(coord_to_id, key)
            coord_to_id[key] = i
            push!(uniq_ids, i)
        end
    end
    nu = length(uniq_ids)
    nu < 3 && return NTuple{3,Int}[]

    # Working point arrays: 1..nu are the distinct input points (in uniq_ids
    # order); nu+1..nu+3 are the super-triangle vertices.
    px = Vector{Float64}(undef, nu + 3)
    py = Vector{Float64}(undef, nu + 3)
    xmn = Inf; xmx = -Inf; ymn = Inf; ymx = -Inf
    @inbounds for q in 1:nu
        i = uniq_ids[q]
        x = X[1, i]; y = X[2, i]
        px[q] = x; py[q] = y
        x < xmn && (xmn = x); x > xmx && (xmx = x)
        y < ymn && (ymn = y); y > ymx && (ymx = y)
    end

    # Super-triangle: large margin around the bbox. Must enclose every point.
    dx = xmx - xmn; dy = ymx - ymn
    dmax = max(dx, dy)
    dmax <= 0 && (dmax = 1.0)            # all points coincident already filtered
    cx = (xmn + xmx) / 2; cy = (ymn + ymx) / 2
    M = 1000 * dmax                      # generous margin
    px[nu + 1] = cx - 2M; py[nu + 1] = cy - M
    px[nu + 2] = cx + 2M; py[nu + 2] = cy - M
    px[nu + 3] = cx;      py[nu + 3] = cy + 2M

    st1 = nu + 1; st2 = nu + 2; st3 = nu + 3

    s = _DelaunayState(px, py,
                       NTuple{3,Int}[], NTuple{3,Int}[], Bool[], Int[])
    sizehint!(s.tri_v, 2 * nu + 8)
    sizehint!(s.tri_adj, 2 * nu + 8)
    sizehint!(s.tri_dead, 2 * nu + 8)

    # Seed with the super-triangle (CCW: st1,st2,st3). All edges are boundary.
    _new_triangle!(s, (st1, st2, st3), (0, 0, 0))

    # Insertion order: Hilbert sort the distinct points (ids 1..nu).
    order = _hilbert_order(px, py, collect(1:nu))

    # Reused scratch buffers for the cavity flood-fill.
    cavity = Int[]; bnd_u = Int[]; bnd_v = Int[]; bnd_out = Int[]
    sizehint!(cavity, 64); sizehint!(bnd_u, 64)
    sizehint!(bnd_v, 64); sizehint!(bnd_out, 64)

    last_t = 1
    @inbounds for p in order
        last_t = _insert_point!(s, p, last_t, cavity, bnd_u, bnd_v, bnd_out)
    end

    # Collect triangles not touching the super-triangle, remap to ORIGINAL X
    # column ids via uniq_ids.
    out = NTuple{3,Int}[]
    sizehint!(out, length(s.tri_v))
    @inbounds for t in 1:length(s.tri_v)
        s.tri_dead[t] && continue
        v = s.tri_v[t]
        (v[1] > nu || v[2] > nu || v[3] > nu) && continue   # touches super-tri
        out_i = uniq_ids[v[1]]; out_j = uniq_ids[v[2]]; out_k = uniq_ids[v[3]]
        push!(out, (out_i, out_j, out_k))
    end
    return out
end
