# DBSCAN backend.
#
# Density-Based Spatial Clustering of Applications with Noise, via Clustering.jl.
# Configuration subtypes `AbstractClusterConfig`; `cluster(smld, ::DBSCANConfig)`
# writes per-emitter cluster labels into `emitter.id` (0 = noise, 1..K = cluster)
# and returns `(smld_out, ClusterInfo)`.

using Clustering

"""
    DBSCANConfig(; eps_nm, min_points=5, use_3d=false, per_dataset=true, remove_unclustered=false)

Configuration for DBSCAN clustering of SMLM localizations.

# Fields
- `eps_nm::Float64`: neighborhood radius in **nanometers**. Coordinates on
  `AbstractEmitter` subtypes are in microns; the backend converts internally.
- `min_points::Int = 5`: minimum number of points in an ε-neighborhood for a
  point to be a core point (classical DBSCAN `minPts`). Also used as the
  minimum cluster size.
- `use_3d::Bool = false`: if `true`, cluster in (x, y, z); requires
  `Emitter3DFit` emitters. Otherwise cluster in (x, y).
- `per_dataset::Bool = true`: if `true`, cluster within each `dataset` index
  independently so that `(dataset, id)` uniquely identifies a cluster in a
  multi-dataset SMLD. If `false`, all emitters are clustered together and
  `id` alone identifies the cluster.
- `remove_unclustered::Bool = false`: if `true`, emitters tagged as noise
  (`id == 0`) are dropped from the returned SMLD.

# Example
```julia
cfg = DBSCANConfig(eps_nm=50.0, min_points=5)
(smld_out, info) = cluster(smld, cfg)
```

See also: [`AbstractClusterConfig`](@ref), [`ClusterInfo`](@ref), [`cluster`](@ref).
"""
Base.@kwdef struct DBSCANConfig <: AbstractClusterConfig
    eps_nm::Float64
    min_points::Int = 5
    use_3d::Bool = false
    per_dataset::Bool = true
    remove_unclustered::Bool = false
end

# Build a d×n matrix of emitter coordinates in microns. `use_3d` requires
# Emitter3DFit (or any emitter that has a :z field); callers get a clear error
# otherwise.
function _coords_matrix(emitters::AbstractVector{<:SMLMData.AbstractEmitter}, use_3d::Bool)
    n = length(emitters)
    if use_3d
        isempty(emitters) || hasproperty(first(emitters), :z) ||
            error("DBSCANConfig(use_3d=true) requires 3D emitters (e.g. Emitter3DFit); " *
                  "got $(eltype(emitters)).")
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

function cluster(smld::SMLMData.BasicSMLD, cfg::DBSCANConfig)
    t0 = time()
    n_in = length(smld.emitters)
    cfg.eps_nm > 0 || throw(ArgumentError("DBSCANConfig.eps_nm must be > 0 (got $(cfg.eps_nm))"))
    cfg.min_points >= 1 ||
        throw(ArgumentError("DBSCANConfig.min_points must be ≥ 1 (got $(cfg.min_points))"))

    # Radius in microns (emitter coordinates are microns; eps is nm).
    radius_μm = cfg.eps_nm / 1000.0

    # Group emitter indices by dataset if per_dataset, otherwise one global group.
    groups = if cfg.per_dataset
        buckets = Dict{Int, Vector{Int}}()
        @inbounds for (i, e) in pairs(smld.emitters)
            push!(get!(() -> Int[], buckets, e.dataset), i)
        end
        # Iterate datasets in sorted order for deterministic label numbering.
        [buckets[k] for k in sort!(collect(keys(buckets)))]
    else
        [collect(1:n_in)]
    end

    cluster_sizes = Int[]
    n_clustered = 0

    for idxs in groups
        # Skip empty groups defensively; Clustering.dbscan on 0 points is fine but
        # we avoid building a zero-column matrix.
        isempty(idxs) && continue

        sub = view(smld.emitters, idxs)
        X = _coords_matrix(sub, cfg.use_3d)
        res = Clustering.dbscan(
            X, radius_μm;
            min_neighbors = cfg.min_points,
            min_cluster_size = cfg.min_points,
        )

        # res.assignments has one entry per column of X, in the order of `idxs`.
        # Labels are 0 (noise) or 1..K within this group; we write them directly
        # to emitter.id, so per-dataset label namespaces stay local per V3.
        @inbounds for (j, i) in pairs(idxs)
            smld.emitters[i].id = res.assignments[j]
        end
        append!(cluster_sizes, res.counts)
        n_clustered += sum(res.counts)
    end

    n_clusters = length(cluster_sizes)
    n_noise = n_in - n_clustered

    # Output SMLD: honor remove_unclustered by filtering emitters with id == 0.
    out_emitters = cfg.remove_unclustered ?
        [e for e in smld.emitters if e.id != 0] :
        smld.emitters
    smld_out = SMLMData.BasicSMLD(
        out_emitters,
        smld.camera,
        smld.n_frames,
        smld.n_datasets,
        smld.metadata,
    )

    info = ClusterInfo(
        n_in,
        n_clustered,
        n_noise,
        n_clusters,
        cluster_sizes,
        :dbscan,
        time() - t0,
    )
    return smld_out, info
end
