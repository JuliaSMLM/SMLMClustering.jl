# Abstract interface and result type for the clustering backends.
# Concrete backends (DBSCAN, HDBSCAN, Voronoi, hierarchical) live in
# their own files and add `cluster` methods dispatched on their own
# concrete subtype of `AbstractClusterConfig`.

"""
    AbstractClusterConfig <: SMLMData.AbstractSMLMConfig

Abstract supertype for SMLMClustering backend configurations.

Each backend defines a concrete subtype (e.g. `DBSCANConfig`,
`HDBSCANConfig`, `VoronoiConfig`, `HierarchicalConfig`) and adds a
`cluster(smld, cfg)` method specialized on it.

# Shared fields (expected on every concrete subtype)
- `min_points::Int`: minimum points for a valid cluster
- `use_3d::Bool`: whether the z-coordinate participates in clustering
- `per_dataset::Bool`: if `true`, cluster within each `dataset` independently
  so that `(dataset, id)` uniquely identifies a cluster across a multi-dataset SMLD
- `remove_unclustered::Bool`: if `true`, drop emitters with `id == 0` from the
  returned SMLD (noise/rejected localizations)

Algorithm-specific fields (`eps_nm`, `min_cluster_size`, `density_factor`,
`cut_nm`, ...) live on the concrete subtype that needs them.
"""
abstract type AbstractClusterConfig <: SMLMData.AbstractSMLMConfig end

"""
    ClusterInfo <: SMLMData.AbstractSMLMInfo

Secondary output from `cluster()` — summary of a clustering run.

# Fields
- `n_locs_in::Int`: number of input localizations
- `n_clustered::Int`: number of localizations assigned to a cluster (`id > 0`)
- `n_noise::Int`: number of localizations tagged as noise (`id == 0`)
- `n_clusters::Int`: number of distinct clusters formed
- `cluster_sizes::Vector{Int}`: size of each cluster, indexed by cluster id
  (`cluster_sizes[k]` is the size of cluster `k`); length equals `n_clusters`
- `algorithm::Symbol`: backend identifier (`:dbscan`, `:hdbscan`, `:voronoi`, `:hierarchical`)
- `elapsed_s::Float64`: wall-clock time spent in the `cluster` call, in seconds

# Example
```julia
(smld_out, info) = cluster(smld, DBSCANConfig(eps_nm=50.0, min_points=5))
println("\$(info.n_clustered)/\$(info.n_locs_in) clustered into \$(info.n_clusters) clusters")
```
"""
struct ClusterInfo <: SMLMData.AbstractSMLMInfo
    n_locs_in::Int
    n_clustered::Int
    n_noise::Int
    n_clusters::Int
    cluster_sizes::Vector{Int}
    algorithm::Symbol
    elapsed_s::Float64
end

"""
    cluster(smld::SMLMData.BasicSMLD, cfg::AbstractClusterConfig) -> (smld, ClusterInfo)

Cluster the localizations in `smld` using the backend selected by the
concrete type of `cfg`. Each backend mutates / rewrites `emitter.id` so
that `0` marks noise and `1..K` mark distinct clusters; if
`cfg.remove_unclustered` is `true` the returned SMLD contains only
clustered emitters.

Dispatch on `AbstractClusterConfig` alone has no implementation —
concrete backends supply a method specialized on their own config type.
Calling `cluster` with an unsupported config raises a clear error.

See also: `AbstractClusterConfig`, `ClusterInfo`.
"""
function cluster(smld::SMLMData.BasicSMLD, cfg::AbstractClusterConfig)
    error("SMLMClustering.cluster has no method for config type $(typeof(cfg)). " *
          "Each backend (DBSCAN, HDBSCAN, Voronoi, hierarchical) adds its own " *
          "`cluster(smld, cfg::SomeConfig)` method; use one of those concrete config types.")
end
