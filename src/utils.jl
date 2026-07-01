# Shared internal helpers — coordinate extraction, distance computation,
# dataset grouping, cluster-label compaction, and output SMLD construction.
# Loaded before backend files; all functions are package-private (underscore prefix).

# Build a d×n coordinate matrix in microns from a vector of emitters.
# d = 2 when use_3d=false; d = 3 when use_3d=true (requires :z property).
function _coords_matrix(emitters::AbstractVector{<:SMLMData.AbstractEmitter}, use_3d::Bool)
    n = length(emitters)
    if use_3d
        isempty(emitters) || hasproperty(first(emitters), :z) ||
            error("use_3d=true requires 3D emitters (e.g. Emitter3DFit); got $(eltype(emitters)).")
        X = Matrix{Float64}(undef, 3, n)
        @inbounds for i in 1:n
            e = emitters[i]
            X[1, i] = e.x
            X[2, i] = e.y
            X[3, i] = e.z
        end
        return X
    else
        X = Matrix{Float64}(undef, 2, n)
        @inbounds for i in 1:n
            e = emitters[i]
            X[1, i] = e.x
            X[2, i] = e.y
        end
        return X
    end
end

# Symmetric n×n pairwise Euclidean distance matrix from a d×n column-major matrix.
_pairwise_distances(X::Matrix{Float64}) = Distances.pairwise(Euclidean(), X; dims=2)

# Group emitter indices by dataset (sorted, deterministic). When
# per_dataset=false, returns a single all-indices group so downstream code
# can always iterate `for idxs in groups`.
function _group_by_dataset(smld::SMLMData.BasicSMLD, per_dataset::Bool)
    n = length(smld.emitters)
    per_dataset || return [collect(1:n)]
    buckets = Dict{Int,Vector{Int}}()
    @inbounds for (i, e) in pairs(smld.emitters)
        push!(get!(() -> Int[], buckets, e.dataset), i)
    end
    [buckets[k] for k in sort!(collect(keys(buckets)))]
end

# Given raw component sizes indexed by 1..k_raw, produce a compact label map
# (raw → final id, or 0 for below-threshold components), push kept sizes onto
# `cluster_sizes`, and return (label_map, n_added). The label_map is local to
# the current group so per-dataset namespaces stay separate (V3).
function _compact_relabel!(cluster_sizes::Vector{Int},
                           raw_counts::Vector{Int},
                           min_points::Int)
    k_raw = length(raw_counts)
    label_map = zeros(Int, k_raw)
    k_local = 0
    added = 0
    @inbounds for (orig, cnt) in enumerate(raw_counts)
        if cnt >= min_points
            k_local += 1
            label_map[orig] = k_local
            push!(cluster_sizes, cnt)
            added += cnt
        end
    end
    label_map, added
end

# Build the output SMLD, honoring `remove_unclustered` by dropping emitters
# with `id == 0`. Camera, frame/dataset counts, and metadata are preserved
# from the input SMLD.
function _build_output(smld::SMLMData.BasicSMLD, remove_unclustered::Bool)
    out_emitters = remove_unclustered ?
        [e for e in smld.emitters if e.id != 0] :
        smld.emitters
    SMLMData.BasicSMLD(out_emitters, smld.camera, smld.n_frames,
                       smld.n_datasets, smld.metadata)
end

# Per-emitter Voronoi cell areas (μm²) for a vector of 2D emitters.
# Returns `(areas, tri)` where `areas` is a Vector{Float64} of length n in the
# same order as `emitters`, and `tri` is the underlying DelaunayTriangulation
# (so callers can reuse it for the Delaunay-adjacency edge set).
#
# Behavior:
# - n < 3: returns `(fill(NaN, n), nothing)` — no tessellation possible.
# - duplicate (x,y): raises `ArgumentError` before triangulation (mirrors
#   `VoronoiConfig`'s guard — duplicate generators cause `get_area` to throw
#   `KeyError`).
# - 3D not supported here; this helper is 2D-only (DelaunayTriangulation.jl
#   limitation, V7). Callers must validate `use_3d == false` upstream.
#
# Cells are clipped to the convex hull (`voronoi(tri; clip=true)`) so every
# generator has a finite area; hull cells are smaller than their infinite-plane
# area (V8 boundary-handling caveat).
function _voronoi_areas(emitters::AbstractVector{<:SMLMData.AbstractEmitter})
    n = length(emitters)
    if n < 3
        return (fill(NaN, n), nothing)
    end
    pts = [(emitters[j].x, emitters[j].y) for j in 1:n]  # μm
    length(unique(pts)) == n ||
        throw(ArgumentError(
            "Voronoi-density helper: group of $n points contains duplicate " *
            "(x,y) coordinates; deduplicate input localizations before " *
            "calling this backend."))
    tri = DelaunayTriangulation.triangulate(pts)
    vor = DelaunayTriangulation.voronoi(tri; clip = true)
    areas = Vector{Float64}(undef, n)
    @inbounds for j in 1:n
        areas[j] = DelaunayTriangulation.get_area(vor, j)
    end
    return (areas, tri)
end

# Even-odd (ray-casting) point-in-polygon test for a closed 2D polygon given as a
# vector of (x, y) vertices. Shared by EdgeClassify (boundary classification) and
# the Hopkins backend (region-restricted reference sampling).
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
