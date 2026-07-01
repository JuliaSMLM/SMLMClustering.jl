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

function Base.show(io::IO, info::ClusterInfo)
    print(io, "ClusterInfo(",
          info.n_clustered, "/", info.n_locs_in, " clustered, ",
          info.n_clusters, " clusters, algorithm=:", info.algorithm, ", ",
          round(info.elapsed_s * 1e3, digits = 1), " ms)")
end

"""
    cluster(smld::SMLMData.BasicSMLD, cfg::AbstractClusterConfig) -> (smld_out, ClusterInfo)

Cluster the localizations in `smld` using the backend selected by the
concrete type of `cfg`. The input `smld` is **not modified** — each backend
deep-copies the input emitters and writes cluster labels onto the copy's
`emitter.id` (`0` marks noise, `1..K` mark distinct clusters). If
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

# ============================================================================
# cluster_statistics: sibling interface for diagnostic / spatial-stat backends
# ============================================================================

"""
    AbstractStatisticsConfig <: SMLMData.AbstractSMLMConfig

Abstract supertype for `cluster_statistics` backend configurations.

Sibling to `AbstractClusterConfig`. Concrete subtypes (e.g. `HopkinsConfig`)
configure spatial-statistic computations that **do not modify** the SMLD —
they read the coordinate set, compute a scalar (and optional vector) summary,
and return it alongside the input SMLD reference.

# Shared fields (expected on every concrete subtype)
- `use_3d::Bool`: whether the z-coordinate participates in the statistic
- `per_dataset::Bool`: if `true`, compute per dataset and aggregate

Algorithm-specific fields (`n_samples`, `seed`, `random_repeats`, ...) live
on the concrete subtype that needs them.
"""
abstract type AbstractStatisticsConfig <: SMLMData.AbstractSMLMConfig end

"""
    ClusterStatisticsInfo <: SMLMData.AbstractSMLMInfo

Secondary output from `cluster_statistics()` — summary of a spatial-statistic
computation.

# Fields
- `n_locs_in::Int`: number of input localizations
- `statistic::Float64`: primary scalar result (Hopkins H, median density, ...)
- `statistic_name::Symbol`: identifier for `statistic` (`:hopkins`, `:median_density`, ...)
- `algorithm::Symbol`: backend identifier (`:hopkins`, `:voronoi_density`, ...)
- `elapsed_s::Float64`: wall-clock time spent in the `cluster_statistics` call (seconds)
- `extras::Dict{Symbol,Any}`: per-backend supplementary outputs — vector-valued
  results (e.g. per-dataset Hopkins scores under `:hopkins_per_dataset`,
  per-emitter densities under `:density_per_emitter`) live here

# Convention for vector-valued backends
Backends that produce a natural vector output (per-dataset Hopkins, per-emitter
density, per-cluster silhouette) place the vector in `extras` under a descriptive
key and a meaningful summary scalar (mean, median, ...) in `statistic`. This keeps
the simple `info.statistic` access ergonomic for one-number consumers while
preserving the full result for callers that need it.

# Example
```julia
(_, info) = cluster_statistics(smld, HopkinsConfig(n_samples=50, seed=1))
println("Hopkins H = \$(round(info.statistic, digits=3))")
per_ds = info.extras[:hopkins_per_dataset]  # Vector{Float64} when per_dataset=true
```
"""
struct ClusterStatisticsInfo <: SMLMData.AbstractSMLMInfo
    n_locs_in::Int
    statistic::Float64
    statistic_name::Symbol
    algorithm::Symbol
    elapsed_s::Float64
    extras::Dict{Symbol,Any}
end

function Base.show(io::IO, info::ClusterStatisticsInfo)
    print(io, "ClusterStatisticsInfo(",
          info.statistic_name, "=", round(info.statistic, digits = 4),
          ", n_locs_in=", info.n_locs_in,
          ", algorithm=:", info.algorithm, ", ",
          round(info.elapsed_s * 1e3, digits = 1), " ms)")
end

"""
    cluster_statistics(smld::SMLMData.BasicSMLD, cfg::AbstractStatisticsConfig)
        -> (smld, ClusterStatisticsInfo)

Compute a spatial / clustering statistic on `smld` using the backend selected
by the concrete type of `cfg`. Returns a tuple `(smld, info)`.

# Pass-through SMLD semantic (NOT non-mutating copy)

Unlike `cluster()` — which deep-copies emitters so it can write cluster labels
without touching the input — `cluster_statistics()` writes nothing onto the
SMLD and returns the **same reference** as the input. The two-element tuple
shape is preserved for ecosystem symmetry (every SMLM step returns
`(smld, info)`), but callers should know the SMLD is the unmodified input,
not a fresh copy. This asymmetry is intentional: copying for a read-only
operation is wasted work, and statistic backends never have a label to
write back.

Dispatch on `AbstractStatisticsConfig` alone has no implementation —
concrete backends supply a method specialized on their own config type.
Calling `cluster_statistics` with an unsupported config raises a clear error
naming the available concrete backends.

See also: `AbstractStatisticsConfig`, `ClusterStatisticsInfo`,
[`cluster`](@ref) (the labeling sibling).
"""
function cluster_statistics(smld::SMLMData.BasicSMLD, cfg::AbstractStatisticsConfig)
    error("SMLMClustering.cluster_statistics has no method for config type $(typeof(cfg)). " *
          "Each backend (e.g. HopkinsConfig, VoronoiDensityConfig) adds its own " *
          "`cluster_statistics(smld, cfg::SomeConfig)` method; use one of those concrete config types.")
end
