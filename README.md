# SMLMClustering

Clustering backends for single-molecule localization microscopy (SMLM) data.

Provides a single entry point over three algorithm families — DBSCAN,
Voronoi-tessellation (SR-Tesseler-style), and agglomerative hierarchical — operating
on `SMLMData.BasicSMLD` datasets.

## Entry point

```julia
(smld_out, info) = cluster(smld, cfg)
```

`cfg` is a concrete config struct that selects the backend and carries its parameters.
The call mutates `emitter.id` on the input SMLD so that `0` marks noise and `1..K`
marks distinct cluster ids. When `cfg.remove_unclustered = true` the returned `smld_out`
contains only clustered emitters; otherwise it shares the input's emitter vector.
`info` is a `ClusterInfo` summary (counts, sizes, algorithm, elapsed time).

## Backends

### DBSCAN

Distance-based density clustering via [Clustering.jl](https://github.com/JuliaStats/Clustering.jl).

```julia
cfg = DBSCANConfig(
    eps_nm          = 50.0,   # neighborhood radius in nm (required)
    min_points      = 5,      # core-point threshold / min cluster size
    use_3d          = false,  # include z-coordinate
    per_dataset     = true,   # cluster within each dataset independently
    remove_unclustered = false,
)
(smld_out, info) = cluster(smld, cfg)
```

Good default choice: scales to large datasets, no O(n²) memory, works in 2D and 3D.

### Voronoi (SR-Tesseler)

Density clustering via Voronoi tessellation, following
[Levet et al., Nat. Methods 2015](https://doi.org/10.1038/nmeth.3579).
A localization is "dense" when its Voronoi cell area is smaller than
`mean_area / density_factor`; dense Delaunay-adjacent points form clusters.

```julia
cfg = VoronoiConfig(
    density_factor  = 2.0,   # density threshold multiplier
    min_points      = 5,     # minimum cluster size
    per_dataset     = true,
    remove_unclustered = false,
)
(smld_out, info) = cluster(smld, cfg)
```

**2D only.** `use_3d = true` raises `ArgumentError`.
Groups with fewer than 3 points are tagged all-noise.
Groups containing exact-duplicate (x,y) coordinates raise `ArgumentError`.

### Hierarchical

Agglomerative hierarchical clustering via `Clustering.hclust` + `cutree`.

```julia
cfg = HierarchicalConfig(
    cut_nm          = 200.0,  # dendrogram cut height in nm (required)
    linkage         = :ward,  # :single | :complete | :average | :ward
    min_points      = 5,
    use_3d          = false,
    per_dataset     = true,
    remove_unclustered = false,
)
(smld_out, info) = cluster(smld, cfg)
```

Builds an O(n²) pairwise distance matrix per group — prefer DBSCAN for
datasets with ≫10,000 localizations per group. Supports 2D and 3D.

**Note on Ward linkage and `cut_nm`:** Ward's dendrogram heights are
variance-increase costs (μm²), not Euclidean distances. The `cut_nm` field
converts to μm before cutting (`h = cut_nm / 1000.0`), but under Ward the
numerical scale is a merging cost, not a distance. Distance-based linkages
(`:single`, `:complete`, `:average`) respect the nanometer semantics directly.

## ClusterInfo fields

| Field | Type | Meaning |
|-------|------|---------|
| `n_locs_in` | `Int` | Input localization count |
| `n_clustered` | `Int` | Localizations assigned to a cluster (`id > 0`) |
| `n_noise` | `Int` | Noise localizations (`id == 0`) |
| `n_clusters` | `Int` | Number of distinct clusters |
| `cluster_sizes` | `Vector{Int}` | Size of each cluster, indexed by cluster id |
| `algorithm` | `Symbol` | `:dbscan`, `:voronoi`, or `:hierarchical` |
| `elapsed_s` | `Float64` | Wall-clock time of the `cluster` call (seconds) |

## Shared config fields

Every backend config struct carries these fields with the same defaults:

| Field | Default | Meaning |
|-------|---------|---------|
| `min_points` | `5` | Minimum points for a valid cluster |
| `use_3d` | `false` | Include z-coordinate in clustering |
| `per_dataset` | `true` | Cluster within each dataset independently |
| `remove_unclustered` | `false` | Drop noise emitters from output |

When `per_dataset = true`, `(dataset, id)` uniquely identifies a cluster across a
multi-dataset SMLD. Cluster ids are local to each dataset; the same id in different
datasets refers to different clusters.

## Installation

```julia
# From the JuliaSMLM GitHub org (once the repo is public):
using Pkg
Pkg.add(url="https://github.com/JuliaSMLM/SMLMClustering.jl")
```

## Dependencies

- [SMLMData.jl](https://github.com/JuliaSMLM/SMLMData.jl) — emitter types and SMLD container
- [Clustering.jl](https://github.com/JuliaStats/Clustering.jl) — DBSCAN and hierarchical clustering
- [Distances.jl](https://github.com/JuliaStats/Distances.jl) — pairwise distance matrix
- [DelaunayTriangulation.jl](https://github.com/JuliaGeometry/DelaunayTriangulation.jl) — Voronoi tessellation
