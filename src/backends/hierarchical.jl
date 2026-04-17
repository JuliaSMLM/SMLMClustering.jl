# Hierarchical clustering backend.
#
# Agglomerative hierarchical clustering via Clustering.hclust + cutree.
# `HierarchicalConfig` subtypes `AbstractClusterConfig`; `cluster(smld, ::HierarchicalConfig)`
# builds an O(n²) pairwise distance matrix, cuts the dendrogram at `cut_nm`, then
# relabels clusters smaller than `min_points` as noise (id = 0).

"""
    HierarchicalConfig(; cut_nm, linkage=:ward, min_points=5, use_3d=false,
                        per_dataset=true, remove_unclustered=false)

Configuration for agglomerative hierarchical clustering of SMLM localizations.

# Fields
- `cut_nm::Float64`: dendrogram cut height in **nanometers** (unit convention V6).
  The dendrogram is cut at `h = cut_nm / 1000.0` μm.
- `linkage::Symbol = :ward`: linkage criterion — `:single`, `:complete`, `:average`,
  or `:ward`. Ward minimizes within-cluster variance and is a reasonable default.
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
cfg = HierarchicalConfig(cut_nm=200.0, linkage=:ward)
(smld_out, info) = cluster(smld, cfg)
```

See also: [`AbstractClusterConfig`](@ref), [`ClusterInfo`](@ref), [`cluster`](@ref).
"""
Base.@kwdef struct HierarchicalConfig <: AbstractClusterConfig
    cut_nm::Float64
    linkage::Symbol = :ward
    min_points::Int = 5
    use_3d::Bool = false
    per_dataset::Bool = true
    remove_unclustered::Bool = false
end

function cluster(smld::SMLMData.BasicSMLD, cfg::HierarchicalConfig)
    t0 = time_ns()
    n_in = length(smld.emitters)
    cfg.cut_nm > 0 ||
        throw(ArgumentError("HierarchicalConfig.cut_nm must be > 0 (got $(cfg.cut_nm))"))
    cfg.min_points >= 1 ||
        throw(ArgumentError("HierarchicalConfig.min_points must be ≥ 1 (got $(cfg.min_points))"))
    cfg.linkage in (:single, :complete, :average, :ward) ||
        throw(ArgumentError(
            "HierarchicalConfig.linkage must be :single, :complete, :average, or :ward " *
            "(got $(cfg.linkage))"))

    cut_h = cfg.cut_nm / 1000.0  # nm → μm
    groups = _group_by_dataset(smld, cfg.per_dataset)

    cluster_sizes = Int[]
    n_clustered = 0

    for idxs in groups
        isempty(idxs) && continue

        sub = view(smld.emitters, idxs)
        X = _coords_matrix(sub, cfg.use_3d)

        D = _pairwise_distances(X)
        hc = Clustering.hclust(D, linkage = cfg.linkage)
        raw_labels = Clustering.cutree(hc, h = cut_h)  # 1..K_raw, no zeros

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
