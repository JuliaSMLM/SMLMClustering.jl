```@meta
CurrentModule = SMLMClustering
```

# SMLMClustering

Clustering and spatial-statistic backends for single-molecule localization
microscopy (SMLM) point data, built on
[`SMLMData.BasicSMLD`](https://github.com/JuliaSMLM/SMLMData.jl).

The package exposes **three parallel verbs**, each dispatched on a concrete
configuration type so that the algorithm is selected by *which* config you pass:

| Verb | Purpose | Returns |
|------|---------|---------|
| [`cluster`](@ref) | assign a cluster label to every emitter | `(smld_out, ::ClusterInfo)` |
| [`cluster_statistics`](@ref) | compute read-only spatial statistics | `(smld, ::ClusterStatisticsInfo)` |
| [`classify_emitters`](@ref) | edge / membrane / interior classification | `(smld, ::EdgeClassifyInfo)` |

See the [User Guide](@ref "User Guide") for the calling conventions and the
[Methods overview](@ref "Methods overview") for the full backend catalog and the
concepts behind each one.

![The same localization field clustered by six backends](assets/comparison_grid.png)

*One synthetic localization field run through six labeling backends — see the [Methods overview](@ref "Methods overview").*

## Quick start

```julia
using SMLMClustering

# Density clustering with DBSCAN: every emitter gets a label on a deep copy.
cfg = DBSCANConfig(eps_nm = 50.0, min_points = 5)
smld_out, info = cluster(smld, cfg)

info.n_clusters        # number of clusters found
info.n_noise           # emitters left unclustered (label 0)
info.cluster_sizes     # size of each cluster, indexed by id

# Read-only spatial statistic: is there clustering tendency at all?
_, stats = cluster_statistics(smld, HopkinsConfig())
stats.statistic        # Hopkins H (≈0.5 random, →1 clustered)

# Edge / membrane / interior classification (2D):
_, edge = classify_emitters(smld, KdeValleyConfig())
edge.class             # Vector{Symbol}: :outside / :membrane / :interior
```

`cluster` is **non-mutating** — input emitters are deep-copied and labels are
written onto the copy's `emitter.id` (`0` = noise, `1..K` = clusters).
`cluster_statistics` is **pass-through** — it returns the *same* SMLD reference
unchanged alongside an info struct. `classify_emitters` returns a **new** SMLD (the
input's metadata copied, with the per-emitter class added under
`"edge_classify_class"`) alongside its info.

## Backend catalog

**Labeling** (`cluster`): [DBSCAN](@ref), [HDBSCAN](@ref), [Hierarchical](@ref),
[Voronoi (SR-Tesseler)](@ref "Voronoi (SR-Tesseler)"),
[MRF density-regime](@ref "MRF density-regime"),
[Point hysteresis](@ref "Point hysteresis").

**Spatial statistics** (`cluster_statistics`): [Hopkins statistic](@ref),
[Voronoi density](@ref), [Local contrast](@ref).

**Edge classification** (`classify_emitters`):
[Edge / Membrane Classification](@ref) (`OuterPolygonConfig`, `KdeValleyConfig`).

## Installation

```julia
using Pkg
Pkg.add("SMLMClustering")          # once registered in the General registry
```

Until then, install directly from GitHub:

```julia
Pkg.add(url = "https://github.com/JuliaSMLM/SMLMClustering.jl")
```

## License

MIT. Developed in the [Lidke Lab](https://github.com/JuliaSMLM) at the University
of New Mexico.
