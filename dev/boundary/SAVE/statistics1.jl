# Cluster geometry statistics derived from boundary polygons.
#
# `boundary_statistics` computes per-cluster area and perimeter from the
# polygon coordinates stored in a `ClusterBoundaryInfo` returned by
# `boundary_clusters`.  It is the first step in a general "cluster statistics
# from coordinates and boundaries" pipeline; future extensions (centroid,
# eccentricity, aspect ratio, ...) will add fields here.
#
# The helpers `_polygon_area` and `_polygon_perimeter` operate on closed
# coordinate vectors directly (no index indirection), consistent with the
# closed polygons stored in `ClusterBoundaryInfo.boundaries_x/y`.

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

# ---------------------------------------------------------------------------
# ClusterBoundaryStatisticsInfo
# ---------------------------------------------------------------------------

"""
    ClusterBoundaryStatisticsInfo

Result returned as the second element of [`boundary_statistics`](@ref).

# Fields
- `n_clusters::Int`: number of clusters for which statistics were computed.
- `cluster_keys::Vector{Tuple{Int,Int}}`: sorted `(dataset, cluster_id)` key
  for each entry — identical ordering to the `ClusterBoundaryInfo` input.
- `areas::Vector{Float64}`: polygon area (μm²) for each cluster.
  `NaN` for clusters whose boundary polygon has fewer than 3 unique vertices.
- `perimeters::Vector{Float64}`: polygon perimeter (μm) for each cluster.
  `NaN` for the same degenerate cases as `areas`.
- `elapsed_s::Float64`: wall-clock time in seconds.

# Accessing statistics
```julia
# By sorted position (1-indexed):
j = 2
println("Area = ", sinfo.areas[j], " μm², Perimeter = ", sinfo.perimeters[j], " μm")

# By (dataset, cluster_id) key:
j = findfirst(==((1, 3)), sinfo.cluster_keys)
println("Area = ", sinfo.areas[j])
```

See also: [`boundary_statistics`](@ref), [`ClusterBoundaryInfo`](@ref).
"""
struct ClusterBoundaryStatisticsInfo
    n_clusters::Int
    cluster_keys::Vector{Tuple{Int,Int}}
    areas::Vector{Float64}
    perimeters::Vector{Float64}
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
    boundary_statistics(smld_out, binfo::ClusterBoundaryInfo)
        -> (smld_out, ClusterBoundaryStatisticsInfo)

Compute geometric statistics — area and perimeter — for every cluster whose
boundary polygon is stored in `binfo`.

# Arguments
- `smld_out`: a labeled `BasicSMLD` (typically the first return value of
  [`cluster`](@ref)).  Accepted for future extension (e.g., per-emitter
  statistics keyed to each cluster); currently unused.
- `binfo`: a [`ClusterBoundaryInfo`](@ref) returned by
  [`boundary_clusters`](@ref).

# Returns
A tuple `(smld_out, info)` where:
- `smld_out` is the **same reference** as the input — read-only, no copy.
- `info` is a [`ClusterBoundaryStatisticsInfo`](@ref) with per-cluster
  `areas` (μm²) and `perimeters` (μm) in the same order as
  `binfo.cluster_keys`.

# Degenerate polygons
Clusters whose boundary polygon has fewer than 3 unique vertices (i.e.
`length(boundaries_x[j]) < 4`) receive `NaN` for both area and perimeter.

# Example
```julia
(smld_out, _) = cluster(smld, DBSCANConfig(eps_nm = 50.0, min_points = 5))
(_, binfo)    = boundary_clusters(smld_out, BoundaryClustersConfig())
(_, sinfo)    = boundary_statistics(smld_out, binfo)

for j in 1:sinfo.n_clusters
    println("Cluster ", sinfo.cluster_keys[j],
            ": area = ", round(sinfo.areas[j], digits = 4), " μm²",
            ", perimeter = ", round(sinfo.perimeters[j], digits = 4), " μm")
end
```

See also: [`ClusterBoundaryStatisticsInfo`](@ref), [`boundary_clusters`](@ref),
[`ClusterBoundaryInfo`](@ref).
"""
function boundary_statistics(smld_out::SMLMData.BasicSMLD,
                             binfo::ClusterBoundaryInfo)
    t0 = time_ns()
    n  = binfo.n_clusters

    areas      = Vector{Float64}(undef, n)
    perimeters = Vector{Float64}(undef, n)

    for j in 1:n
        bx = binfo.boundaries_x[j]
        by = binfo.boundaries_y[j]
        if length(bx) < 4   # fewer than 3 unique vertices
            areas[j]      = NaN
            perimeters[j] = NaN
        else
            areas[j]      = _polygon_area(bx, by)
            perimeters[j] = _polygon_perimeter(bx, by)
        end
    end

    info = ClusterBoundaryStatisticsInfo(
        n,
        binfo.cluster_keys,
        areas,
        perimeters,
        (time_ns() - t0) / 1e9,
    )

    return smld_out, info
end
