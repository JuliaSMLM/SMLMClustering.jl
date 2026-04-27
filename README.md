# SMLMClustering

Clustering and spatial-statistic backends for single-molecule localization
microscopy (SMLM) data.

Provides two parallel entry points over `SMLMData.BasicSMLD`:

- **`cluster`** — labeling backends (DBSCAN, Voronoi-tessellation /
  SR-Tesseler-style, agglomerative hierarchical) that assign cluster ids to
  every emitter.
- **`cluster_statistics`** — read-only spatial-statistic backends (Hopkins
  clustering tendency, Voronoi per-emitter density, ...) that compute
  summary scalars / vectors.

## Entry points

```julia
# Labeling
(smld_out, info) = cluster(smld, cfg)

# Read-only spatial statistics
(smld_passthrough, stats_info) = cluster_statistics(smld, stats_cfg)
```

`cluster()` is **non-mutating**: input emitters are deep-copied, cluster labels
are written onto the copy's `emitter.id` (`0` marks noise, `1..K` mark distinct
clusters), and `info::ClusterInfo` carries the size summary. When
`cfg.remove_unclustered = true` the returned `smld_out` contains only clustered
emitters.

`cluster_statistics()` is **pass-through**: it returns the same SMLD reference
as input (no allocation, no mutation) alongside a
`ClusterStatisticsInfo` summary. The two-tuple shape is preserved for ecosystem
symmetry, but callers should treat the first element as the unmodified input.

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
# Distance-based linkage: cut_threshold is in nm.
cfg = HierarchicalConfig(
    cut_threshold   = 200.0,   # cut height; unit depends on linkage (see below)
    linkage         = :single, # :single | :complete | :average | :ward
    min_points      = 5,
    use_3d          = false,
    per_dataset     = true,
    remove_unclustered = false,
)
(smld_out, info) = cluster(smld, cfg)

# Ward linkage: specify number of clusters directly (units-agnostic).
cfg_ward = HierarchicalConfig(n_clusters = 3, linkage = :ward)
```

Exactly one of `cut_threshold` or `n_clusters` must be supplied; providing both
or neither raises `ArgumentError`.

Builds an O(n²) pairwise distance matrix per group — prefer DBSCAN for
datasets with ≫10,000 localizations per group. Supports 2D and 3D.

**Unit convention for `cut_threshold`:** for distance-based linkages
(`:single`, `:complete`, `:average`) the value is in **nanometers** and is
converted to μm internally. For `:ward` the dendrogram height is a
variance-increase cost (roughly μm²) and is passed through without conversion —
there is no meaningful nm interpretation under Ward, which is why `n_clusters`
is usually the cleaner choice for Ward.

## Spatial-statistic backend (`cluster_statistics`)

### Hopkins

Hopkins clustering-tendency statistic, sample-based with a KDTree NN backend
(`NearestNeighbors.jl`).

```julia
cfg = HopkinsConfig(
    n_samples       = 50,        # reference / sampled point count per repeat
    random_repeats  = 10,        # average over independent repeats
    seed            = 1,         # RNG seed for reproducibility
    use_3d          = false,
    per_dataset     = true,      # report per-dataset H in extras + mean as `statistic`
)
(_, info) = cluster_statistics(smld, cfg)
println("Hopkins H = ", round(info.statistic, digits = 3))
println("per-dataset H = ", info.extras[:hopkins_per_dataset])
```

**Interpretation:**
- `H ≈ 0.5`: data consistent with uniform spatial randomness (Poisson)
- `H → 1.0`: strong clustering tendency
- `H → 0.0`: anti-clustering / regular spacing

Edge cases (empty group, `n_samples > n_points`, zero-extent bbox) return
`NaN` for the affected group rather than erroring.

### Voronoi density

Per-emitter Voronoi cell area and corresponding local density `ρᵢ = 1/Aᵢ`,
intended for downstream thresholding (Otsu / GMM on `log ρ`, fixed cutoff,
etc.) — e.g. cell-structure masking on dense membrane regions.

```julia
cfg = VoronoiDensityConfig(
    use_3d      = false,    # 2D only — same constraint as VoronoiConfig
    per_dataset = true,     # tessellate each dataset independently
)
(_, info) = cluster_statistics(smld, cfg)
println("median density = ", round(info.statistic, digits = 2), " μm⁻²")

ρ = info.extras[:density_per_emitter]   # Vector{Float64}, length == n_locs_in
A = info.extras[:area_per_emitter]      # Vector{Float64}, μm²
# ρ[i] and A[i] correspond to smld.emitters[i] (flat in original emitter order,
# NOT grouped by dataset).
```

**2D only.** `use_3d = true` raises `ArgumentError`.
Groups with fewer than 3 points: those emitters get `NaN` density / area.
Groups containing exact-duplicate (x,y) coordinates raise `ArgumentError`.

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

## ClusterStatisticsInfo fields

| Field | Type | Meaning |
|-------|------|---------|
| `n_locs_in` | `Int` | Input localization count |
| `statistic` | `Float64` | Primary scalar result (e.g. Hopkins H) |
| `statistic_name` | `Symbol` | Identifier for `statistic` (`:hopkins`, ...) |
| `algorithm` | `Symbol` | Backend identifier (`:hopkins`, ...) |
| `elapsed_s` | `Float64` | Wall-clock time of the `cluster_statistics` call (seconds) |
| `extras` | `Dict{Symbol,Any}` | Per-backend supplementary outputs (vector results, per-group breakdowns, ...) |

**Convention for vector-valued backends:** put a meaningful summary scalar in
`statistic` (mean, median, ...) and the full vector in `extras` under a
descriptive key (e.g. Hopkins per-dataset vector under
`:hopkins_per_dataset`). This keeps the simple `info.statistic` access
ergonomic while preserving the full result.

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
- [NearestNeighbors.jl](https://github.com/KristofferC/NearestNeighbors.jl) — KDTree NN queries (Hopkins backend)
