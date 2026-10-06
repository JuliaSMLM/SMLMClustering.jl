# Cluster geometry and spacing statistics derived from boundary polygons and
# raw member-point coordinates.
#
# `boundary_statistics` computes, for every cluster in a `ClusterBoundaryInfo`
# returned by `boundary_clusters`, the full statistic set defined in
# `dev/boundary/STATISTICS`: area/perimeter at three boundary tightnesses
# (the `binfo.shrink` boundary, the convex hull, and the tightest alpha
# shape), shape descriptors derived from those (equivalent radius,
# compactness, circularity, convexity, solidity), point-count statistics,
# and inter-point / inter-cluster spacing statistics (sigma_actual,
# cluster_width, nearest-neighbor spacing, center-to-center and
# edge-to-edge distances between clusters).
#
# `n_clusters`, `n_points`, `n_clustered`, and `n_isolated` from
# `dev/boundary/STATISTICS` are intentionally not duplicated here — they are
# already available from `ClusterInfo` (`n_clusters`, `n_locs_in`,
# `n_clustered`, `n_noise`), the secondary output of `cluster()`.
#
# The helpers `_polygon_area` and `_polygon_perimeter` operate on closed
# coordinate vectors directly (no index indirection), consistent with the
# closed polygons stored in `ClusterBoundaryInfo.boundaries_x/y`.
#
# Reference:
# Wooten ZT, Yu C, Court LE, Peterson CB. "Predictive modeling using shape
# statistics for interpretable and robust quality assurance of automated
# contours in radiation treatment planning." PMCID: PMC10091357,
# NIHMSID: NIHMS1888738, PMID: 36540994.
# https://pmc.ncbi.nlm.nih.gov/articles/PMC10091357/

# ---------------------------------------------------------------------------
# Internal geometry helpers
# ---------------------------------------------------------------------------

# Shoelace formula on closed coordinate vectors (bx[1] == bx[end]).
# Returns the unsigned area — always ≥ 0.
function _polygon_area(bx::AbstractVector{<:Real}, by::AbstractVector{<:Real})
    n = length(bx) - 1   # last element repeats first
    A = 0.0
    @inbounds for k in 1:n
        A += bx[k] * by[k + 1] - bx[k + 1] * by[k]
    end
    return abs(A) * 0.5
end

# Sum of edge lengths around the closed polygon.
function _polygon_perimeter(bx::AbstractVector{<:Real}, by::AbstractVector{<:Real})
    n = length(bx) - 1
    P = 0.0
    @inbounds for k in 1:n
        dx = bx[k + 1] - bx[k]
        dy = by[k + 1] - by[k]
        P += sqrt(dx^2 + dy^2)
    end
    return P
end

# Area/perimeter of the polygon traced by index vector `idx` (closed,
# idx[1] == idx[end]) into coordinate vectors `xc`/`yc`. NaN for fewer than
# 3 unique vertices, matching the existing binfo-boundary degenerate rule.
function _indexed_area_perimeter(xc::AbstractVector{Float64}, yc::AbstractVector{Float64},
                                  idx::Vector{Int})
    length(idx) < 4 && return (NaN, NaN)
    bx = @view xc[idx]
    by = @view yc[idx]
    return (_polygon_area(bx, by), _polygon_perimeter(bx, by))
end

# All pairwise Euclidean distances among n ≥ 2 points.
function _pairwise_distances(x::AbstractVector{Float64}, y::AbstractVector{Float64})
    n = length(x)
    d = Vector{Float64}(undef, n * (n - 1) ÷ 2)
    k = 0
    @inbounds for i in 1:(n - 1), j in (i + 1):n
        k += 1
        d[k] = hypot(x[i] - x[j], y[i] - y[j])
    end
    return d
end

# Maximum pairwise distance among the unique points in (x, y); NaN if < 2.
function _max_pairwise_dist(x::AbstractVector{Float64}, y::AbstractVector{Float64})
    n = length(x)
    n < 2 && return NaN
    dmax = 0.0
    @inbounds for i in 1:(n - 1), j in (i + 1):n
        d = hypot(x[i] - x[j], y[i] - y[j])
        d > dmax && (dmax = d)
    end
    return dmax
end

# Minimum distance over all point pairs drawn from two point sets; NaN if
# either set is empty. Used as the vertex-pair approximation to the
# edge-to-edge distance between two clusters' convex hulls.
function _min_dist_between_point_sets(x1::AbstractVector{Float64}, y1::AbstractVector{Float64},
                                       x2::AbstractVector{Float64}, y2::AbstractVector{Float64})
    (isempty(x1) || isempty(x2)) && return NaN
    dmin = Inf
    @inbounds for i in eachindex(x1), j in eachindex(x2)
        d = hypot(x1[i] - x2[j], y1[i] - y2[j])
        d < dmin && (dmin = d)
    end
    return dmin
end

# Nearest same-cluster-neighbor distance for every point in (x, y).
# n < 2 ⇒ no such distance exists ⇒ returns an empty vector.
function _nn_within(x::Vector{Float64}, y::Vector{Float64})
    n = length(x)
    n < 2 && return Float64[]
    X = permutedims(hcat(x, y))                # 2 × n, required by KDTree
    tree = NearestNeighbors.KDTree(X)
    nn = Vector{Float64}(undef, n)
    @inbounds for i in 1:n
        _, dists = NearestNeighbors.knn(tree, view(X, :, i), 2, true)
        nn[i] = dists[2]                        # dists[1] == 0 (self)
    end
    return nn
end

# Group emitter indices by cluster key, exactly mirroring `boundary_clusters`'
# own grouping rule so that lookups by `binfo.cluster_keys` land on the same
# raw points that produced `binfo`'s boundary polygons.
function _group_emitter_indices(emitters, per_dataset::Bool)
    groups = Dict{Tuple{Int,Int}, Vector{Int}}()
    @inbounds for (i, e) in pairs(emitters)
        e.id == 0 && continue
        key = per_dataset ? (e.dataset, e.id) : (1, e.id)
        push!(get!(Vector{Int}, groups, key), i)
    end
    return groups
end

# ---------------------------------------------------------------------------
# ClusterBoundaryStatisticsInfo
# ---------------------------------------------------------------------------

"""
    ClusterBoundaryStatisticsInfo

Result returned as the second element of [`boundary_statistics`](@ref).

Every field below is a `Vector` in the same order as `cluster_keys` (one
entry per cluster), except `min_c2c_dist`, `min_e2e_dist`, and `elapsed_s`,
which are single overall scalars.

# Identity
- `n_clusters::Int`: number of clusters for which statistics were computed.
- `cluster_keys::Vector{Tuple{Int,Int}}`: sorted `(dataset, cluster_id)` key
  for each entry — identical ordering to the `ClusterBoundaryInfo` input.

# Boundary indices (local — index into that cluster's own point array, not
  the SMLD)
- `indices::Vector{Vector{Int}}`: closed boundary polygon at `binfo.shrink`.
- `indicesConvex::Vector{Vector{Int}}`: closed convex-hull boundary (shrink=0).
- `indicesTight::Vector{Vector{Int}}`: closed tightest alpha-shape boundary
  (shrink=1).
All three are `Int[]` for clusters with no matching points in `smld_out`
(see "Missing raw points" below).

# Area / perimeter (μm² / μm)
- `area`, `perimeter`: at `binfo.shrink`, computed directly from
  `binfo.boundaries_x/y` (unchanged from the original `boundary_statistics`).
- `areaConvex`, `perimeterConvex`: convex-hull boundary.
- `areaTight`, `perimeterTight`: tightest alpha-shape boundary.
NaN wherever the underlying polygon has fewer than 3 unique vertices.

# Shape descriptors (dimensionless)
- `equiv_radius = sqrt(area / π) / 2`.
- `compactness = 4π·area / perimeter²` (1 for a circle, 0 for a line segment).
- `circularity = 4π·areaTight / perimeterConvex²` (1 for a circle; relatively
  insensitive to irregular boundaries).
- `convexity = perimeterConvex / perimeterTight` (1 for a convex object; < 1
  if the boundary is irregular).
- `solidity = areaTight / areaConvex` (1 for a solid object; < 1 if the
  boundary is irregular or has holes).

# Point counts and spacing
- `center::Vector{Tuple{Float64,Float64}}`: (x, y) centroid of each cluster's
  member points, μm.
- `n_pts_per_cluster::Vector{Int}`: number of member points per cluster.
- `n_pts_per_area::Vector{Float64}`: `n_pts_per_cluster / area`; NaN for
  clusters with fewer than 3 points.
- `sigma_actual::Vector{Float64}`: standard deviation of all pairwise
  distances between member points within the cluster; NaN for < 2 points.
- `cluster_width::Vector{Float64}`: maximum pairwise distance among the
  cluster's boundary vertices (`indices`); NaN for < 2 unique vertices.
- `nn_within_clusters::Vector{Vector{Float64}}`: nearest same-cluster-neighbor
  distance for each member point (length == `n_pts_per_cluster[j]`); empty
  for clusters with < 2 points.

# Inter-cluster distances (μm)
- `min_c2c_dists::Vector{Float64}`: for each cluster, the minimum
  center-to-center distance to any other cluster.
- `min_e2e_dists::Vector{Float64}`: for each cluster, the minimum
  edge-to-edge distance to any other cluster, approximated as the nearest
  pair of convex-hull vertices between the two clusters.
- `min_c2c_dist::Float64`, `min_e2e_dist::Float64`: `minimum` of the above
  two vectors over all valid entries. NaN when fewer than 2 clusters have
  valid centers/hulls.

- `elapsed_s::Float64`: wall-clock time in seconds.

# Missing raw points
`area`/`perimeter`/`equiv_radius`/`compactness` need only `binfo`'s boundary
polygon and are always computed. Every other field needs the cluster's raw
member points, looked up in `smld_out` by `(dataset, id)`; if a cluster key
has no matching emitters (e.g. `smld_out` is unrelated to `binfo`, as in
some tests), those fields are `NaN` / empty for that cluster rather than
raising an error.

# Accessing statistics
```julia
# By sorted position (1-indexed):
j = 2
println("Area = ", sinfo.area[j], " μm², compactness = ", sinfo.compactness[j])

# By (dataset, cluster_id) key:
j = findfirst(==((1, 3)), sinfo.cluster_keys)
println("Area = ", sinfo.area[j])
```

See also: [`boundary_statistics`](@ref), [`ClusterBoundaryInfo`](@ref).
"""
struct ClusterBoundaryStatisticsInfo
    n_clusters::Int
    cluster_keys::Vector{Tuple{Int,Int}}
    area::Vector{Float64}
    perimeter::Vector{Float64}
    indices::Vector{Vector{Int}}
    indicesConvex::Vector{Vector{Int}}
    indicesTight::Vector{Vector{Int}}
    areaConvex::Vector{Float64}
    areaTight::Vector{Float64}
    perimeterConvex::Vector{Float64}
    perimeterTight::Vector{Float64}
    center::Vector{Tuple{Float64,Float64}}
    equiv_radius::Vector{Float64}
    compactness::Vector{Float64}
    circularity::Vector{Float64}
    convexity::Vector{Float64}
    solidity::Vector{Float64}
    n_pts_per_cluster::Vector{Int}
    n_pts_per_area::Vector{Float64}
    sigma_actual::Vector{Float64}
    cluster_width::Vector{Float64}
    min_c2c_dists::Vector{Float64}
    min_e2e_dists::Vector{Float64}
    min_c2c_dist::Float64
    min_e2e_dist::Float64
    nn_within_clusters::Vector{Vector{Float64}}
    elapsed_s::Float64
end

function Base.show(io::IO, info::ClusterBoundaryStatisticsInfo)
    print(io, "ClusterBoundaryStatisticsInfo(",
          info.n_clusters, " clusters, ",
          round(info.elapsed_s * 1e3, digits = 1), " ms)")
end

# ---------------------------------------------------------------------------
# boundary_statistics: public API
# ---------------------------------------------------------------------------

"""
    boundary_statistics(smld_out, binfo::ClusterBoundaryInfo; per_dataset=true)
        -> (smld_out, ClusterBoundaryStatisticsInfo)

Compute the full geometry and spacing statistic set (see
[`ClusterBoundaryStatisticsInfo`](@ref)) for every cluster whose boundary
polygon is stored in `binfo`.

# Arguments
- `smld_out`: a labeled `BasicSMLD` (typically the first return value of
  [`cluster`](@ref)), used to look up each cluster's raw member-point
  coordinates by `(dataset, id)`.
- `binfo`: a [`ClusterBoundaryInfo`](@ref) returned by
  [`boundary_clusters`](@ref).
- `per_dataset`: must match the value used in the `BoundaryClustersConfig`
  that produced `binfo` (and, before that, in the `cluster()` call), so that
  raw points are regrouped identically to how `binfo`'s cluster keys were
  built.

# Returns
A tuple `(smld_out, info)` where:
- `smld_out` is the **same reference** as the input — read-only, no copy.
- `info` is a [`ClusterBoundaryStatisticsInfo`](@ref).

# Example
```julia
(smld_out, _) = cluster(smld, DBSCANConfig(eps_nm = 50.0, min_points = 5))
(_, binfo)    = boundary_clusters(smld_out, BoundaryClustersConfig())
(_, sinfo)    = boundary_statistics(smld_out, binfo)

for j in 1:sinfo.n_clusters
    println("Cluster ", sinfo.cluster_keys[j],
            ": area = ", round(sinfo.area[j], digits = 4), " μm²",
            ", compactness = ", round(sinfo.compactness[j], digits = 3))
end
```

See also: [`ClusterBoundaryStatisticsInfo`](@ref), [`boundary_clusters`](@ref),
[`ClusterBoundaryInfo`](@ref).
"""
function boundary_statistics(smld_out::SMLMData.BasicSMLD,
                              binfo::ClusterBoundaryInfo;
                              per_dataset::Bool = true)
    t0 = time_ns()
    n  = binfo.n_clusters

    area      = Vector{Float64}(undef, n)
    perimeter = Vector{Float64}(undef, n)

    indices       = Vector{Vector{Int}}(undef, n)
    indicesConvex = Vector{Vector{Int}}(undef, n)
    indicesTight  = Vector{Vector{Int}}(undef, n)

    areaConvex      = Vector{Float64}(undef, n)
    areaTight       = Vector{Float64}(undef, n)
    perimeterConvex = Vector{Float64}(undef, n)
    perimeterTight  = Vector{Float64}(undef, n)

    center            = Vector{Tuple{Float64,Float64}}(undef, n)
    n_pts_per_cluster = Vector{Int}(undef, n)
    n_pts_per_area    = Vector{Float64}(undef, n)
    sigma_actual      = Vector{Float64}(undef, n)
    cluster_width     = Vector{Float64}(undef, n)
    nn_within_clusters   = Vector{Vector{Float64}}(undef, n)

    # Convex-hull coordinates per cluster, retained for the second pass that
    # computes inter-cluster edge-to-edge distances.
    hull_x = Vector{Vector{Float64}}(undef, n)
    hull_y = Vector{Vector{Float64}}(undef, n)

    groups = _group_emitter_indices(smld_out.emitters, per_dataset)

    for j in 1:n
        # --- area/perimeter at binfo.shrink: unchanged, from binfo's own polygon ---
        bx = binfo.boundaries_x[j]
        by = binfo.boundaries_y[j]
        if length(bx) < 4
            area[j]      = NaN
            perimeter[j] = NaN
        else
            area[j]      = _polygon_area(bx, by)
            perimeter[j] = _polygon_perimeter(bx, by)
        end

        # --- raw member points for this cluster (may be empty) ---
        idxs = get(groups, binfo.cluster_keys[j], Int[])
        npts = length(idxs)
        n_pts_per_cluster[j] = npts

        if npts == 0
            indices[j] = indicesConvex[j] = indicesTight[j] = Int[]
            areaConvex[j] = areaTight[j] = perimeterConvex[j] = perimeterTight[j] = NaN
            center[j]          = (NaN, NaN)
            n_pts_per_area[j]  = NaN
            sigma_actual[j]    = NaN
            cluster_width[j]   = NaN
            nn_within_clusters[j] = Float64[]
            hull_x[j] = hull_y[j] = Float64[]
        else
            emitters = smld_out.emitters
            xc = Vector{Float64}(undef, npts)
            yc = Vector{Float64}(undef, npts)
            @inbounds for (k, i) in pairs(idxs)
                xc[k] = emitters[i].x
                yc[k] = emitters[i].y
            end

            indices[j]       = boundary_cluster(xc, yc; shrink = binfo.shrink)
            indicesConvex[j] = boundary_cluster(xc, yc; shrink = 0.0)
            indicesTight[j]  = boundary_cluster(xc, yc; shrink = 1.0)

            areaConvex[j], perimeterConvex[j] = _indexed_area_perimeter(xc, yc, indicesConvex[j])
            areaTight[j],  perimeterTight[j]  = _indexed_area_perimeter(xc, yc, indicesTight[j])

            center[j]         = (mean(xc), mean(yc))
            n_pts_per_area[j] = npts >= 3 ? npts / area[j] : NaN
            sigma_actual[j]   = npts >= 2 ? std(_pairwise_distances(xc, yc)) : NaN

            verts = unique(indices[j])
            cluster_width[j] = _max_pairwise_dist(xc[verts], yc[verts])

            nn_within_clusters[j] = _nn_within(xc, yc)

            hull_x[j] = xc[indicesConvex[j]]
            hull_y[j] = yc[indicesConvex[j]]
        end
    end

    # --- shape descriptors: derived from area/perimeter, always computed ---
    equiv_radius = @. sqrt(area / pi) / 2
    compactness = @. 4 * pi * area / perimeter^2
    circularity = @. 4 * pi * areaTight / perimeterConvex^2
    convexity   = @. perimeterConvex / perimeterTight
    solidity    = @. areaTight / areaConvex

    # --- inter-cluster distances: second pass, needs every center/hull first ---
    min_c2c_dists = fill(NaN, n)
    min_e2e_dists = fill(NaN, n)
    for j in 1:n
        cxj, cyj = center[j]
        isnan(cxj) && continue
        best_c2c = Inf
        best_e2e = Inf
        for jp in 1:n
            jp == j && continue
            cxjp, cyjp = center[jp]
            if !isnan(cxjp)
                d = hypot(cxj - cxjp, cyj - cyjp)
                d < best_c2c && (best_c2c = d)
            end
            if !isempty(hull_x[j]) && !isempty(hull_x[jp])
                d = _min_dist_between_point_sets(hull_x[j], hull_y[j], hull_x[jp], hull_y[jp])
                d < best_e2e && (best_e2e = d)
            end
        end
        min_c2c_dists[j] = isfinite(best_c2c) ? best_c2c : NaN
        min_e2e_dists[j] = isfinite(best_e2e) ? best_e2e : NaN
    end
    valid_c2c    = filter(!isnan, min_c2c_dists)
    valid_e2e    = filter(!isnan, min_e2e_dists)
    min_c2c_dist = isempty(valid_c2c) ? NaN : minimum(valid_c2c)
    min_e2e_dist = isempty(valid_e2e) ? NaN : minimum(valid_e2e)

    info = ClusterBoundaryStatisticsInfo(
        n,
        binfo.cluster_keys,
        area,
        perimeter,
        indices,
        indicesConvex,
        indicesTight,
        areaConvex,
        areaTight,
        perimeterConvex,
        perimeterTight,
        center,
        equiv_radius,
        compactness,
        circularity,
        convexity,
        solidity,
        n_pts_per_cluster,
        n_pts_per_area,
        sigma_actual,
        cluster_width,
        min_c2c_dists,
        min_e2e_dists,
        min_c2c_dist,
        min_e2e_dist,
        nn_within_clusters,
        (time_ns() - t0) / 1e9,
    )

    return smld_out, info
end
