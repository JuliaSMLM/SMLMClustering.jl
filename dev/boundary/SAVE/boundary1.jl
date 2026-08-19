# Cluster boundary computation: alpha-shape boundary extraction for 2D SMLM data.
#
# `boundary_cluster` computes the boundary polygon of a single 2D point cloud,
# analogous to MATLAB's `boundary(x, y, shrink)`.  The algorithm builds a
# Delaunay triangulation (via DelaunayTriangulation.jl, already a project
# dependency), then selects a subset of triangles whose circumradius falls
# below a data-adaptive threshold controlled by `shrink`:
#
#   shrink = 0  →  alpha = max circumradius  →  all triangles included
#                  ⇒ result equals the convex hull
#   shrink = 1  →  alpha = min circumradius  →  only the smallest triangles
#                  ⇒ tightest possible boundary
#
# Boundary edges are edges that belong to exactly one included triangle
# (interior edges belong to two; hull edges of excluded regions belong to
# zero).  The boundary polygon is traced from these edges and returned as a
# closed, counter-clockwise `Vector{Int}` of indices into the input x/y
# arrays, with the first index repeated at the end — the same convention
# as MATLAB's `boundary`.
#
# `boundary_clusters` applies `boundary_cluster` to every cluster in a
# labeled SMLD and returns per-cluster boundary polygons as (x, y) coordinate
# vectors in microns.
#
# Reference: Edelsbrunner, Kirkpatrick, Seidel, "On the shape of a set of
# points in the plane," IEEE Trans. Inform. Theory 29(4), 1983.

# ---------------------------------------------------------------------------
# Internal geometry helpers
# ---------------------------------------------------------------------------

# Circumradius of the triangle with vertices (ax,ay), (bx,by), (cx,cy).
# Uses R = abc / (2 |cross|) = abc / (4 Area).
# Returns Inf for degenerate (zero-area) triangles so they sort to the end
# and are automatically excluded from any finite alpha complex.
@inline function _circumradius_2d(ax::Float64, ay::Float64,
                                   bx::Float64, by::Float64,
                                   cx::Float64, cy::Float64)
    # Side lengths
    ab = sqrt((bx - ax)^2 + (by - ay)^2)
    bc = sqrt((cx - bx)^2 + (cy - by)^2)
    ca = sqrt((ax - cx)^2 + (ay - cy)^2)
    # |2 × Area| via the z-component of the cross product
    twice_area = abs((bx - ax) * (cy - ay) - (cx - ax) * (by - ay))
    # Guard: treat near-zero-area triangles as degenerate.
    # Threshold scaled by perimeter² keeps the check dimensionally consistent.
    twice_area < 1e-14 * (ab * bc + bc * ca + ca * ab + 1e-30) && return Inf
    return (ab * bc * ca) / (2.0 * twice_area)
end

# Signed area of a closed polygon given by an index sequence into x/y.
# Positive ⇒ counter-clockwise; negative ⇒ clockwise (standard math convention).
# The loop vector must satisfy loop[1] == loop[end] (closed polygon).
function _signed_area_polygon(x::AbstractVector{<:Real},
                               y::AbstractVector{<:Real},
                               loop::Vector{Int})
    A = 0.0
    n = length(loop) - 1   # last element repeats first
    @inbounds for k in 1:n
        i = loop[k]
        j = loop[k + 1]
        A += x[i] * y[j] - x[j] * y[i]
    end
    return A * 0.5
end

# ---------------------------------------------------------------------------
# Polygon tracing from a boundary-edge adjacency list
# ---------------------------------------------------------------------------

# Trace a single closed loop starting at `start`.  Visited vertices are
# removed from `unvisited` in place.  Each vertex in a valid boundary has
# exactly two neighbors; `prev` is used to avoid immediately back-tracking.
function _trace_single_loop!(adj::Dict{Int,Vector{Int}},
                              start::Int,
                              unvisited::Set{Int})
    loop = Int[start]
    delete!(unvisited, start)

    prev = -1
    v    = start

    # Upper bound prevents infinite loops from malformed adjacency.
    for _ in 1:(length(adj) + 1)
        nbrs = get(adj, v, Int[])
        next = -1

        for w in nbrs
            w == prev && continue           # don't backtrack
            if w == start && length(loop) >= 3
                push!(loop, start)          # close the polygon
                return loop
            end
            if w in unvisited
                next = w
                break
            end
        end

        next == -1 && break                 # dead end — shouldn't happen for valid boundary

        push!(loop, next)
        delete!(unvisited, next)
        prev = v
        v    = next
    end

    # Close the loop even if tracing stopped early (degenerate case).
    isempty(loop) || loop[1] != loop[end] && push!(loop, start)
    return loop
end

# Trace all connected closed loops from the adjacency list.
# Returns a Vector of loops; each loop is a closed index sequence.
function _trace_boundary_loops(adj::Dict{Int,Vector{Int}})
    unvisited = Set{Int}(keys(adj))
    loops     = Vector{Vector{Int}}()

    while !isempty(unvisited)
        start = minimum(unvisited)          # deterministic for reproducibility
        loop  = _trace_single_loop!(adj, start, unvisited)
        length(loop) >= 4 && push!(loops, loop)   # need ≥ 3 unique vertices + closing repeat
    end

    return loops
end

# ---------------------------------------------------------------------------
# boundary_cluster: public low-level API
# ---------------------------------------------------------------------------

"""
    boundary_cluster(x, y; shrink=0.5) -> Vector{Int}

Compute the boundary polygon of a 2D point cloud, analogous to MATLAB's
`boundary(x, y, shrink)`.

# Arguments
- `x`, `y`: coordinate vectors of equal length (in any consistent units).
- `shrink`: shrink factor in [0, 1] controlling boundary tightness.
  - `shrink = 0` → convex hull (loosest possible boundary).
  - `shrink = 1` → alpha shape with the minimum circumradius threshold
    (tightest boundary; only the most compact Delaunay triangles are kept).
  - Intermediate values interpolate between these extremes via the quantile
    of triangle circumradii.

# Returns
A `Vector{Int}` of indices into `x`/`y` tracing the outer boundary polygon
in counter-clockwise order.  The polygon is **closed**: `result[1] == result[end]`.
This matches MATLAB's `boundary` output convention.

# Edge cases
- `n < 3` returns `1:n` (all points; no triangulation possible).
- `n == 3` returns the triangle itself `[1, 2, 3, 1]`.
- Collinear or near-collinear inputs are handled by DelaunayTriangulation.jl.
- Degenerate (zero-area) triangles have circumradius `Inf` and are
  automatically excluded from all non-trivial alpha complexes.

# Algorithm
1. Build the Delaunay triangulation of `(x, y)`.
2. For each solid (non-ghost) triangle, compute its circumradius.
3. Set `alpha_radius = quantile(circumradii, 1 − shrink)`:
   - the `(1−shrink)` quantile of all circumradii.
4. Include triangles with circumradius ≤ `alpha_radius` in the alpha complex.
5. Extract boundary edges (edges in exactly one included triangle).
6. Trace boundary loops; return the largest loop with CCW orientation.

# Performance
O(n log n) dominated by the Delaunay triangulation.  Circumradius computation
and edge counting are O(n) with concrete types, `@inbounds`, and a pre-sized
`Dict`.  Allocation count scales as O(n).

# Example
```julia
using SMLMClustering
x = randn(500); y = randn(500)
k = boundary_cluster(x, y; shrink = 0.5)
# k[i] is an index into x/y; x[k], y[k] traces the boundary polygon.
```

See also: [`boundary_clusters`](@ref), [`BoundaryClustersConfig`](@ref).
"""
function boundary_cluster(x::AbstractVector{<:Real},
                           y::AbstractVector{<:Real};
                           shrink::Float64 = 0.5)
    n = length(x)
    n == length(y) ||
        throw(ArgumentError(
            "boundary_cluster: x and y must have the same length " *
            "(got $(length(x)) and $(length(y)))."))
    0.0 <= shrink <= 1.0 ||
        throw(ArgumentError(
            "boundary_cluster: shrink must be in [0, 1] (got $shrink)."))

    # --- Degenerate small cases ---
    n == 0 && return Int[]
    n == 1 && return [1]
    n == 2 && return [1, 2, 1]
    n == 3 && return [1, 2, 3, 1]

    # Convert to Float64 tuples required by DelaunayTriangulation.
    pts = [(Float64(x[i]), Float64(y[i])) for i in 1:n]
    tri = DelaunayTriangulation.triangulate(pts)

    # --- Convex hull shortcut (shrink = 0) ---
    # get_convex_hull_vertices returns a closed CCW index sequence.
    if shrink == 0.0
        return collect(DelaunayTriangulation.get_convex_hull_vertices(tri))
    end

    # --- Collect solid triangles ---
    # each_solid_triangle yields NTuple{3,Int} with positive vertex indices.
    triangles = NTuple{3,Int}[]
    for t in DelaunayTriangulation.each_solid_triangle(tri)
        push!(triangles, t)
    end
    m = length(triangles)
    m == 0 && return collect(DelaunayTriangulation.get_convex_hull_vertices(tri))

    # --- Circumradius for every solid triangle ---
    circumradii = Vector{Float64}(undef, m)
    @inbounds for (k, (i, j, l)) in pairs(triangles)
        circumradii[k] = _circumradius_2d(
            Float64(x[i]), Float64(y[i]),
            Float64(x[j]), Float64(y[j]),
            Float64(x[l]), Float64(y[l]),
        )
    end

    # --- Alpha radius from shrink factor ---
    # Quantile mapping: shrink=0 → R_max (all triangles) → convex hull.
    #                   shrink=1 → R_min (smallest triangle) → tightest boundary.
    # We sort and index: q_idx = ceil((1−shrink)·m), clamped to [1, m].
    R_sorted = sort(circumradii)
    q_idx    = clamp(ceil(Int, (1.0 - shrink) * m), 1, m)
    alpha_radius = R_sorted[q_idx]

    # --- Count how many in-alpha triangles each canonical edge belongs to ---
    # Boundary edges are those with count == 1.
    # Hull edges of the Delaunay triangulation appear in only 1 solid triangle,
    # so they are automatically counted correctly.
    edge_count = Dict{Tuple{Int,Int}, Int}()
    sizehint!(edge_count, 3 * m)

    @inbounds for (k, (i, j, l)) in pairs(triangles)
        circumradii[k] <= alpha_radius || continue
        # Canonical (min, max) key ensures each edge is counted consistently
        # regardless of the triangle's vertex winding order.
        for (a, b) in ((i, j), (j, l), (l, i))
            e = a < b ? (a, b) : (b, a)
            edge_count[e] = get(edge_count, e, 0) + 1
        end
    end

    # --- Build adjacency list from boundary edges ---
    adj = Dict{Int, Vector{Int}}()
    for ((a, b), cnt) in edge_count
        cnt == 1 || continue
        push!(get!(Vector{Int}, adj, a), b)
        push!(get!(Vector{Int}, adj, b), a)
    end

    isempty(adj) && return collect(DelaunayTriangulation.get_convex_hull_vertices(tri))

    # --- Trace all boundary loops; keep the largest ---
    loops = _trace_boundary_loops(adj)
    isempty(loops) && return collect(DelaunayTriangulation.get_convex_hull_vertices(tri))

    # "Largest" by number of unique boundary vertices (length − 1 because the
    # first vertex is repeated at the end to close the polygon).
    best = loops[1]
    for loop in loops
        length(loop) > length(best) && (best = loop)
    end

    # --- Enforce counter-clockwise orientation ---
    # Positive signed area ⇒ CCW; negative ⇒ CW.
    if _signed_area_polygon(x, y, best) < 0.0 && length(best) > 2
        # Reverse the interior (keep first == last): [v1, vk, ..., v2, v1]
        best = [best[1]; reverse(best[2:end-1]); best[end]]
    end

    return best
end

# ---------------------------------------------------------------------------
# BoundaryClustersConfig and ClusterBoundaryInfo
# ---------------------------------------------------------------------------

"""
    BoundaryClustersConfig(; shrink=0.5, per_dataset=true)

Configuration for [`boundary_clusters`](@ref).

# Fields
- `shrink::Float64 = 0.5`: shrink factor passed to [`boundary_cluster`](@ref)
  for every cluster.  `0` gives a convex-hull boundary; `1` gives the tightest
  alpha-shape boundary.
- `per_dataset::Bool = true`: when `true`, cluster ids in `smld_out` are
  interpreted as dataset-local (i.e., `(dataset, id)` is the unique cluster
  key).  Set to the same value used in the original `cluster()` call.

# Example
```julia
cfg_clust = DBSCANConfig(eps_nm = 50.0, min_points = 5)
(smld_out, _)   = cluster(smld, cfg_clust)

cfg_bnd = BoundaryClustersConfig(shrink = 0.5, per_dataset = true)
(_, binfo) = boundary_clusters(smld_out, cfg_bnd)
# binfo.boundaries_x[j], binfo.boundaries_y[j] are the coordinates of
# the boundary polygon for cluster binfo.cluster_keys[j].
```

See also: [`boundary_clusters`](@ref), [`boundary_cluster`](@ref),
[`ClusterBoundaryInfo`](@ref).
"""
Base.@kwdef struct BoundaryClustersConfig
    shrink::Float64 = 0.5
    per_dataset::Bool = true
end

"""
    ClusterBoundaryInfo

Result returned as the second element of [`boundary_clusters`](@ref).

# Fields
- `n_clusters::Int`: number of clusters for which a boundary was computed.
- `cluster_keys::Vector{Tuple{Int,Int}}`: sorted `(dataset, cluster_id)` key
  for each entry.  When `per_dataset = false` the `dataset` value reflects the
  emitter's actual `.dataset` field, which may be constant or irrelevant for
  indexing.  Use `findfirst(==(key), binfo.cluster_keys)` to locate a specific
  cluster.
- `boundaries_x::Vector{Vector{Float64}}`: x-coordinates (μm) of the boundary
  polygon for each cluster, in the same order as `cluster_keys`.  The polygon
  is closed (`boundaries_x[j][1] == boundaries_x[j][end]`) and CCW.
- `boundaries_y::Vector{Vector{Float64}}`: corresponding y-coordinates (μm).
- `shrink::Float64`: the shrink factor used.
- `elapsed_s::Float64`: wall-clock time in seconds.

# Accessing a cluster's boundary
```julia
# By sorted position (1-indexed):
j = 3
plot(binfo.boundaries_x[j], binfo.boundaries_y[j])

# By (dataset, cluster_id) key:
j = findfirst(==(( 1, 2)), binfo.cluster_keys)
plot(binfo.boundaries_x[j], binfo.boundaries_y[j])
```

See also: [`boundary_clusters`](@ref), [`BoundaryClustersConfig`](@ref).
"""
struct ClusterBoundaryInfo
    n_clusters::Int
    cluster_keys::Vector{Tuple{Int,Int}}
    boundaries_x::Vector{Vector{Float64}}
    boundaries_y::Vector{Vector{Float64}}
    shrink::Float64
    elapsed_s::Float64
end

function Base.show(io::IO, info::ClusterBoundaryInfo)
    print(io, "ClusterBoundaryInfo(",
          info.n_clusters, " clusters, shrink=", info.shrink, ", ",
          round(info.elapsed_s * 1e3, digits = 1), " ms)")
end

# ---------------------------------------------------------------------------
# boundary_clusters: public high-level API
# ---------------------------------------------------------------------------

"""
    boundary_clusters(smld_out, cfg::BoundaryClustersConfig)
        -> (smld_out, ClusterBoundaryInfo)

Compute the boundary polygon for every cluster in a labeled SMLD, using the
alpha-shape algorithm in [`boundary_cluster`](@ref).

# Arguments
- `smld_out`: a `BasicSMLD` whose emitters carry cluster labels in their
  `.id` field (`0` = noise, `1..K` = cluster).  Typically the first return
  value of [`cluster`](@ref).
- `cfg`: a [`BoundaryClustersConfig`](@ref) specifying the shrink factor and
  whether cluster ids are dataset-local.

# Returns
A tuple `(smld_out, info)` where:
- `smld_out` is the **same reference** as the input — this function is
  read-only and does not copy or modify the SMLD (analogous to
  `cluster_statistics`).
- `info` is a [`ClusterBoundaryInfo`](@ref) with per-cluster boundary
  polygons in microns.

# Per-dataset semantics
When `cfg.per_dataset = true`, the unique cluster identifier is the pair
`(emitter.dataset, emitter.id)`.  Set this to match the value used in the
original `cluster()` call so that clusters are not accidentally merged across
datasets.

# Example
```julia
(smld_out, _)   = cluster(smld, DBSCANConfig(eps_nm = 50.0, min_points = 10))
cfg_bnd         = BoundaryClustersConfig(shrink = 0.5)
(_, binfo)      = boundary_clusters(smld_out, cfg_bnd)

for j in 1:binfo.n_clusters
    println("Cluster \$(binfo.cluster_keys[j]): \$(length(binfo.boundaries_x[j])-1) boundary points")
end
```

See also: [`BoundaryClustersConfig`](@ref), [`ClusterBoundaryInfo`](@ref),
[`boundary_cluster`](@ref).
"""
function boundary_clusters(smld_out::SMLMData.BasicSMLD,
                            cfg::BoundaryClustersConfig)
    t0 = time_ns()

    0.0 <= cfg.shrink <= 1.0 ||
        throw(ArgumentError(
            "BoundaryClustersConfig.shrink must be in [0, 1] (got $(cfg.shrink))."))

    emitters = smld_out.emitters

    # --- Group emitter indices by (dataset, cluster_id), skip noise (id == 0) ---
    # When per_dataset=false we still use the emitter's actual .dataset value as
    # the key's first component; cluster ids are then globally unique across
    # datasets (as assigned by cluster()), so no merging occurs.
    groups = Dict{Tuple{Int,Int}, Vector{Int}}()
    @inbounds for (i, e) in pairs(emitters)
        e.id == 0 && continue
        key = cfg.per_dataset ? (e.dataset, e.id) : (1, e.id)
        push!(get!(Vector{Int}, groups, key), i)
    end

    # Sort keys for deterministic, human-readable indexing.
    sorted_keys = sort!(collect(keys(groups)))
    n_clusters  = length(sorted_keys)

    boundaries_x = Vector{Vector{Float64}}(undef, n_clusters)
    boundaries_y = Vector{Vector{Float64}}(undef, n_clusters)

    for (j, key) in pairs(sorted_keys)
        idxs = groups[key]
        n    = length(idxs)

        # Extract 2D coordinates in μm for this cluster.
        xc = Vector{Float64}(undef, n)
        yc = Vector{Float64}(undef, n)
        @inbounds for (k, i) in pairs(idxs)
            xc[k] = emitters[i].x
            yc[k] = emitters[i].y
        end

        # Compute alpha-shape boundary; get back indices into xc/yc.
        bnd_idx = boundary_cluster(xc, yc; shrink = cfg.shrink)

        # Map boundary indices to absolute (x, y) coordinates.
        bx = Vector{Float64}(undef, length(bnd_idx))
        by = Vector{Float64}(undef, length(bnd_idx))
        @inbounds for (k, bi) in pairs(bnd_idx)
            bx[k] = xc[bi]
            by[k] = yc[bi]
        end

        boundaries_x[j] = bx
        boundaries_y[j] = by
    end

    info = ClusterBoundaryInfo(
        n_clusters,
        sorted_keys,
        boundaries_x,
        boundaries_y,
        cfg.shrink,
        (time_ns() - t0) / 1e9,
    )

    return smld_out, info
end
