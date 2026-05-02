# Point-graph hysteresis seed-and-grow clustering.
#
# Connected-component labeling on a kNN graph, gated by:
#   1. seed presence (component must contain ≥1 seed point), and
#   2. minimum component size.
#
# This is the cheaper, more interpretable cousin of the GMM+Potts MRF pipeline
# in `MRFDensityClusterConfig`. The MRF couples unary density evidence with a
# spatial prior implicitly through ICM iterations; hysteresis makes the
# coupling explicit and discrete: high-confidence "seed" points pull in
# adjacent "support" points through a graph BFS, and components without seed
# evidence are dropped wholesale. When the underlying density unary is poorly
# identified (e.g. weak GMM separation, the case that flooded RGY in the
# A431/MAPN test), the MRF amplifies the bad emission model; hysteresis is
# robust because a bad seed mask just means fewer (or no) seeds — not a flood.
#
# Caller supplies the seed and support boolean vectors. Typically these come
# from thresholding `LocalContrastFeature` and/or `_knn_density` outputs; the
# config does not bundle the feature step, so the same backend supports any
# discriminator.

import NearestNeighbors

"""
    point_hysteresis_clusters(smld, seed, support; graph_k=12, min_points=100,
                              use_3d=false, per_dataset=false,
                              remove_unclustered=false)
        -> (smld_out, ClusterInfo)

Cluster `smld` by hysteresis seed-and-grow on a per-emitter kNN graph.

# Algorithm
1. Build a `KDTree` over the (per-dataset or pooled) emitter coordinates.
2. Iterate emitters in input order. Whenever an unvisited support point is
   found, BFS through the support set following kNN edges (`graph_k`
   neighbors per node, excluding self). Mark every reached node visited.
3. After the BFS, if the component contains at least one seed point AND has
   at least `min_points` members, assign it a fresh cluster id (local within
   the per-dataset group, matching the package convention from V3). Otherwise
   leave its members at `id = 0`.

# Arguments
- `smld::BasicSMLD`: input localizations. Not modified — backend deep-copies.
- `seed::AbstractVector{Bool}` of length `length(smld.emitters)`: high-confidence
  foreground points. Components must contain at least one to be kept.
- `support::AbstractVector{Bool}` of length `length(smld.emitters)`: candidate
  foreground points. BFS only crosses through support points.
- `graph_k::Int = 12`: degree of the kNN graph.
- `min_points::Int = 100`: minimum component size for a cluster to be kept.
  Components with `<min_points` members or no seed presence are dropped to
  noise (`id = 0`).
- `use_3d::Bool = false`: build the kNN graph in (x, y, z).
- `per_dataset::Bool = false`: when `true`, cluster within each dataset
  independently. Cluster ids are local to each dataset's namespace so
  `(dataset, id)` is the unique identifier (V3).
- `remove_unclustered::Bool = false`: if `true`, drop emitters with `id = 0`
  from the returned SMLD.

# Convention
- `seed` must imply `support` (every seed point is also a support point).
  Violations raise `ArgumentError`.
- BFS does not cross outside the support set, so isolated seed points without
  support neighbors form components of size 1 and are dropped by `min_points`.

# Example
```julia
# Local-contrast feature → seed/support thresholds → hysteresis
(_, info_f) = cluster_statistics(smld, LocalContrastFeature(density_k=200,
                                                            background_k=2000))
contrast = info_f.extras[:contrast_per_emitter]
fine = info_f.extras[:log_density_per_emitter]
fine_floor = quantile(filter(isfinite, fine), 0.35)
seed = isfinite.(contrast) .& (contrast .> 0.25) .& (fine .> fine_floor)
support = isfinite.(contrast) .& (contrast .> -0.05) .& (fine .> fine_floor)

(smld_out, info) = point_hysteresis_clusters(smld, seed, support;
                                             graph_k=12, min_points=150)
```

See also: [`LocalContrastFeature`](@ref), [`MRFDensityClusterConfig`](@ref)
(the GMM+Potts alternative), [`ClusterInfo`](@ref).
"""
function point_hysteresis_clusters(smld::SMLMData.BasicSMLD,
                                   seed::AbstractVector{Bool},
                                   support::AbstractVector{Bool};
                                   graph_k::Int = 12,
                                   min_points::Int = 100,
                                   use_3d::Bool = false,
                                   per_dataset::Bool = false,
                                   remove_unclustered::Bool = false)
    t0 = time_ns()
    n_in = length(smld.emitters)
    length(seed) == n_in ||
        throw(ArgumentError("seed length $(length(seed)) must equal emitter count $n_in"))
    length(support) == n_in ||
        throw(ArgumentError("support length $(length(support)) must equal emitter count $n_in"))
    graph_k >= 1 ||
        throw(ArgumentError("point_hysteresis_clusters: graph_k must be ≥ 1 (got $graph_k)"))
    min_points >= 1 ||
        throw(ArgumentError("point_hysteresis_clusters: min_points must be ≥ 1 (got $min_points)"))

    @inbounds for i in 1:n_in
        if seed[i] && !support[i]
            throw(ArgumentError(
                "point_hysteresis_clusters: seed implies support, but emitter " *
                "$i is seed=true with support=false. Compute support as a " *
                "superset of seed (e.g. via a looser threshold on the same feature)."))
        end
    end

    # Non-mutating semantics: deep-copy emitters so cluster labels go to a
    # fresh SMLD (V9). Zero all ids on the copy so unclustered points come
    # out as id=0 (the "noise" contract from types.jl) regardless of what
    # the input carried — pipelines that feed the SMLD through prior
    # labeling stages (e.g. BaGoL group ids) would otherwise see those ids
    # leak through on emitters this backend does not classify.
    smld = SMLMData.BasicSMLD(deepcopy(smld.emitters), smld.camera,
                              smld.n_frames, smld.n_datasets, smld.metadata)
    @inbounds for e in smld.emitters
        e.id = 0
    end

    cluster_sizes = Int[]
    n_clustered = 0

    groups = _group_by_dataset(smld, per_dataset)

    for idxs in groups
        n = length(idxs)
        n == 0 && continue

        sub = view(smld.emitters, idxs)
        X = _coords_matrix(sub, use_3d)
        k_use = min(graph_k, n - 1)
        # Skip degenerate group (single point can't form a cluster meeting min_points).
        if k_use < 1 || n < min_points
            continue
        end
        tree = NearestNeighbors.KDTree(X)

        local_seed = falses(n)
        local_support = falses(n)
        @inbounds for (j, i) in pairs(idxs)
            local_seed[j] = seed[i]
            local_support[j] = support[i]
        end

        visited = falses(n)
        stack = Int[]
        comp = Int[]
        n_local_clusters = 0

        @inbounds for start in 1:n
            (local_support[start] && !visited[start]) || continue
            empty!(stack)
            empty!(comp)
            has_seed = false
            visited[start] = true
            push!(stack, start)

            while !isempty(stack)
                v = pop!(stack)
                push!(comp, v)
                local_seed[v] && (has_seed = true)
                nbrs, _ = NearestNeighbors.knn(tree, view(X, :, v), k_use + 1, true)
                for w in nbrs
                    w == v && continue
                    if local_support[w] && !visited[w]
                        visited[w] = true
                        push!(stack, w)
                    end
                end
            end

            if has_seed && length(comp) >= min_points
                n_local_clusters += 1
                push!(cluster_sizes, length(comp))
                for j in comp
                    smld.emitters[idxs[j]].id = n_local_clusters
                    n_clustered += 1
                end
            end
        end
    end

    n_clusters = length(cluster_sizes)
    n_noise = n_in - n_clustered
    smld_out = _build_output(smld, remove_unclustered)

    info = ClusterInfo(
        n_in,
        n_clustered,
        n_noise,
        n_clusters,
        cluster_sizes,
        :point_hysteresis,
        (time_ns() - t0) / 1e9,
    )
    return smld_out, info
end
