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
    t0 = time()
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

    # Group emitter indices by dataset if per_dataset, else one global group.
    groups = if cfg.per_dataset
        buckets = Dict{Int, Vector{Int}}()
        @inbounds for (i, e) in pairs(smld.emitters)
            push!(get!(() -> Int[], buckets, e.dataset), i)
        end
        [buckets[k] for k in sort!(collect(keys(buckets)))]
    else
        [collect(1:n_in)]
    end

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

        # Build remapping: raw cluster → final label, local to this group (V3 namespace).
        # Clusters smaller than min_points become noise (0); the rest get compact 1..K_local.
        label_map = zeros(Int, k_raw)
        k_local = 0
        @inbounds for (orig, cnt) in enumerate(raw_counts)
            if cnt >= cfg.min_points
                k_local += 1
                label_map[orig] = k_local
                push!(cluster_sizes, cnt)
                n_clustered += cnt
            end
        end

        # Write final labels to emitter.id (local namespace per V3).
        @inbounds for (j, i) in pairs(idxs)
            smld.emitters[i].id = label_map[raw_labels[j]]
        end
    end

    n_clusters = length(cluster_sizes)
    n_noise = n_in - n_clustered

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
        :hierarchical,
        time() - t0,
    )
    return smld_out, info
end
