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
