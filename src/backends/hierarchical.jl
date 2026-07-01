# Hierarchical clustering backend.
#
# Agglomerative hierarchical clustering via Clustering.hclust + cutree.
# `HierarchicalConfig` subtypes `AbstractClusterConfig`; `cluster(smld, ::HierarchicalConfig)`
# builds an O(n²) pairwise distance matrix, cuts the dendrogram at `cut_threshold`
# (distance-based linkages) OR at `n_clusters` clusters, then relabels clusters
# smaller than `min_points` as noise (id = 0).

"""
    HierarchicalConfig(; cut_threshold=nothing, n_clusters=nothing, linkage=:ward,
                        min_points=5, use_3d=false, per_dataset=true,
                        remove_unclustered=false)

Configuration for agglomerative hierarchical clustering of SMLM localizations.

Exactly one of `cut_threshold` or `n_clusters` must be supplied.

# Fields
- `cut_threshold::Union{Float64,Nothing} = nothing`: dendrogram cut height.
  **Units depend on linkage:** for distance-based linkages (`:single`, `:complete`,
  `:average`) the value is in **nanometers** and is converted to μm internally
  (`h = cut_threshold / 1000.0`). For `:ward` the dendrogram height is a
  variance-increase cost (roughly μm²) and is passed through without conversion —
  there is no meaningful "nm" interpretation under Ward.
- `n_clusters::Union{Int,Nothing} = nothing`: cut the dendrogram to produce
  exactly this many clusters (before `min_points` filtering). Natural for Ward,
  where `cut_threshold` has no intuitive unit. Mutually exclusive with
  `cut_threshold`.
- `linkage::Symbol = :ward`: linkage criterion — `:single`, `:complete`, `:average`,
  or `:ward`. Ward minimizes within-cluster variance.
- `min_points::Int = 5`: clusters with fewer than `min_points` members after cutting
  are relabeled noise (`id = 0`). Remaining clusters are renumbered compactly `1..K`.
- `use_3d::Bool = false`: include z-coordinate in distance calculation.
- `per_dataset::Bool = true`: cluster within each dataset independently so that
  `(dataset, id)` uniquely identifies a cluster in a multi-dataset SMLD.
- `remove_unclustered::Bool = false`: drop noise emitters (`id == 0`) from the output.

!!! note "Scalability"
    Hierarchical clustering builds an O(n²) pairwise distance matrix. For large datasets
    (≫10,000 localizations per group) DBSCAN is preferred.

# Example
```julia
# Distance-based: cut at 200 nm under single linkage.
cfg = HierarchicalConfig(cut_threshold=200.0, linkage=:single)
(smld_out, info) = cluster(smld, cfg)

# Ward linkage: specify number of clusters directly.
cfg2 = HierarchicalConfig(n_clusters=3, linkage=:ward)
(smld_out, info) = cluster(smld, cfg2)
```

See also: [`AbstractClusterConfig`](@ref), [`ClusterInfo`](@ref), [`cluster`](@ref).
"""
Base.@kwdef struct HierarchicalConfig <: AbstractClusterConfig
    cut_threshold::Union{Float64,Nothing} = nothing
    n_clusters::Union{Int,Nothing} = nothing
    linkage::Symbol = :ward
    min_points::Int = 5
    use_3d::Bool = false
    per_dataset::Bool = true
    remove_unclustered::Bool = false
end

function cluster(smld::SMLMData.BasicSMLD, cfg::HierarchicalConfig)
    t0 = time_ns()
    # Non-mutating semantics: deep-copy emitters so cluster labels go to a
    # fresh SMLD, not back onto the caller's input. See KB V9.
    smld = SMLMData.BasicSMLD(deepcopy(smld.emitters), smld.camera,
                              smld.n_frames, smld.n_datasets, smld.metadata)
    n_in = length(smld.emitters)

    ct_set = cfg.cut_threshold !== nothing
    nc_set = cfg.n_clusters !== nothing
    xor(ct_set, nc_set) ||
        throw(ArgumentError(
            "HierarchicalConfig: exactly one of cut_threshold or n_clusters must be set " *
            "(got cut_threshold=$(cfg.cut_threshold), n_clusters=$(cfg.n_clusters))"))
    !ct_set || cfg.cut_threshold > 0 ||
        throw(ArgumentError(
            "HierarchicalConfig.cut_threshold must be > 0 (got $(cfg.cut_threshold))"))
    !nc_set || cfg.n_clusters >= 1 ||
        throw(ArgumentError(
            "HierarchicalConfig.n_clusters must be ≥ 1 (got $(cfg.n_clusters))"))
    cfg.min_points >= 1 ||
        throw(ArgumentError("HierarchicalConfig.min_points must be ≥ 1 (got $(cfg.min_points))"))
    cfg.linkage in (:single, :complete, :average, :ward) ||
        throw(ArgumentError(
            "HierarchicalConfig.linkage must be :single, :complete, :average, or :ward " *
            "(got $(cfg.linkage))"))

    # Distance-based linkages: cut_threshold is in nm, convert to μm.
    # Ward: cut_threshold is a variance-cost unit (~μm²), pass through.
    cut_h = if ct_set
        cfg.linkage === :ward ? cfg.cut_threshold : cfg.cut_threshold / 1000.0
    else
        nothing
    end

    groups = _group_by_dataset(smld, cfg.per_dataset)

    cluster_sizes = Int[]
    n_clustered = 0

    for idxs in groups
        isempty(idxs) && continue

        sub = view(smld.emitters, idxs)
        X = _coords_matrix(sub, cfg.use_3d)

        D = _pairwise_distances(X)
        hc = Clustering.hclust(D, linkage = cfg.linkage)
        raw_labels = ct_set ?
            Clustering.cutree(hc, h = cut_h) :
            Clustering.cutree(hc, k = min(cfg.n_clusters, length(idxs)))

        k_raw = maximum(raw_labels; init = 0)
        raw_counts = zeros(Int, k_raw)
        @inbounds for l in raw_labels
            raw_counts[l] += 1
        end

        label_map, added = _compact_relabel!(cluster_sizes, raw_counts, cfg.min_points)
        n_clustered += added

        # Write final labels to emitter.id (local namespace per V3).
        @inbounds for (j, i) in pairs(idxs)
            smld.emitters[i].id = label_map[raw_labels[j]]
        end
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
        :hierarchical,
        (time_ns() - t0) / 1e9,
    )
    return smld_out, info
end
