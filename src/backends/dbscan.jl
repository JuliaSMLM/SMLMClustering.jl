# DBSCAN backend.
#
# Density-Based Spatial Clustering of Applications with Noise (Ester et al. 1996).
# Configuration subtypes `AbstractClusterConfig`; `cluster(smld, ::DBSCANConfig)`
# writes per-emitter cluster labels into `emitter.id` (0 = noise, 1..K = cluster)
# and returns `(smld_out, ClusterInfo)`.
#
# Runs on the package's own neighbor graph (`build_precision_neighbor_graph`, every
# pair with d ≤ eps) and core/border pass (`_core_border_labels!`), shared with
# PrecisionDBSCANConfig. It does not wrap `Clustering.dbscan`: Clustering.jl 0.15
# adds the neighbors of every explored point to the cluster before testing whether
# that point is core, so a cluster reached up to 2·eps past its last core point.
#
# Shared helpers (_coords_matrix, _group_by_dataset, etc.) are provided by utils.jl,
# which is included before the backend files.

"""
    DBSCANConfig(; eps_nm, min_points=5, use_3d=false, per_dataset=true, remove_unclustered=false)

Configuration for DBSCAN clustering of SMLM localizations.

Classical DBSCAN (Ester et al. 1996): a point is *core* if at least `min_points`
points, itself included, lie within `eps_nm` of it (distance `≤ eps_nm`). Core points
within `eps_nm` of each other share a cluster. A non-core point joins a cluster only
if it lies within `eps_nm` of one of that cluster's core points (the lowest cluster
id when it is near cores of several clusters); otherwise it is noise. Labels are
deterministic for a given emitter order.

# Fields
- `eps_nm::Float64`: neighborhood radius in **nanometers**. Coordinates on
  `AbstractEmitter` subtypes are in microns; the backend converts internally.
- `min_points::Int = 5`: minimum number of points in an ε-neighborhood for a
  point to be a core point (classical DBSCAN `minPts`). Also used as the
  minimum cluster size: a cluster left with fewer than `min_points` members
  (possible only when neighboring clusters take some of its border points)
  becomes noise.
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
    # Non-mutating semantics: deep-copy emitters so cluster labels go to a
    # fresh SMLD, not back onto the caller's input. See KB V9.
    smld = SMLMData.BasicSMLD(deepcopy(smld.emitters), smld.camera,
                              smld.n_frames, smld.n_datasets, smld.metadata)
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
        isempty(idxs) && continue

        sub = view(smld.emitters, idxs)
        X = _coords_matrix(sub, cfg.use_3d)
        # Every pair with d ≤ eps (KDTree inrange, inclusive), all of them active.
        g = build_precision_neighbor_graph(X, radius_μm)
        labels = Vector{Int}(undef, g.n)
        _core_border_labels!(labels, g, cfg.min_points, _all_active)

        # Each point has one label, so sizes come straight from the labels. Drop
        # clusters below min_points (one can shrink only by losing border points to
        # a lower-id cluster), renumber the survivors 1..k_local and push their
        # sizes onto cluster_sizes. label_map is local per V3.
        group_sizes = zeros(Int, maximum(labels; init = 0))
        @inbounds for a in labels
            a == 0 && continue
            group_sizes[a] += 1
        end
        label_map, added = _compact_relabel!(cluster_sizes, group_sizes, cfg.min_points)
        @inbounds for (j, i) in pairs(idxs)
            a = labels[j]
            smld.emitters[i].id = a == 0 ? 0 : label_map[a]
        end
        n_clustered += added
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
