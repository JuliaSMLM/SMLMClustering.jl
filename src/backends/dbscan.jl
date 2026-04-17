# DBSCAN backend.
#
# Density-Based Spatial Clustering of Applications with Noise, via Clustering.jl.
# Configuration subtypes `AbstractClusterConfig`; `cluster(smld, ::DBSCANConfig)`
# writes per-emitter cluster labels into `emitter.id` (0 = noise, 1..K = cluster)
# and returns `(smld_out, ClusterInfo)`.
#
# Shared helpers (_coords_matrix, _group_by_dataset, etc.) are provided by utils.jl,
# which is included before the backend files.

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

function cluster(smld::SMLMData.BasicSMLD, cfg::DBSCANConfig)
    t0 = time_ns()
    n_in = length(smld.emitters)
    cfg.eps_nm > 0 || throw(ArgumentError("DBSCANConfig.eps_nm must be > 0 (got $(cfg.eps_nm))"))
    cfg.min_points >= 1 ||
        throw(ArgumentError("DBSCANConfig.min_points must be ≥ 1 (got $(cfg.min_points))"))

    # Radius in microns (emitter coordinates are microns; eps is nm).
    radius_μm = cfg.eps_nm / 1000.0

    groups = _group_by_dataset(smld, cfg.per_dataset)

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
    smld_out = _build_output(smld, cfg.remove_unclustered)

    info = ClusterInfo(
        n_in,
        n_clustered,
        n_noise,
        n_clusters,
        cluster_sizes,
        :dbscan,
        (time_ns() - t0) / 1e9,
    )
    return smld_out, info
end
