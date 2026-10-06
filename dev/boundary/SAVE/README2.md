# SMLMClustering

Clustering and spatial-statistic backends for single-molecule localization
microscopy (SMLM) data.

Provides three parallel entry points over `SMLMData.BasicSMLD`:

- **`cluster`** — labeling backends (DBSCAN, Voronoi-tessellation /
  SR-Tesseler-style, agglomerative hierarchical) that assign cluster ids to
  every emitter.
- **`cluster_statistics`** — read-only spatial-statistic backends (Hopkins
  clustering tendency, Voronoi per-emitter density, ...) that compute
  summary scalars / vectors.
- **`classify_emitters`** — edge / membrane / interior classification
  (`OuterPolygonConfig`, and the validated adaptive `KdeValleyConfig`) labeling
  each emitter `:outside` / `:membrane` / `:interior`.

Two further utilities build on `cluster()`'s output rather than dispatching on
their own abstract config family: **`boundary_clusters`** extracts each
cluster's alpha-shape boundary polygon, and **`boundary_statistics`** computes
a full set of per-cluster shape and spacing statistics from that boundary.
**`plot_clusters`** renders a Gaussian super-resolution image with optional
cluster-centroid and boundary-polygon overlays for quick visual QC.

## Entry points

```julia
# Labeling
(smld_out, info) = cluster(smld, cfg)

# Read-only spatial statistics
(smld_passthrough, stats_info) = cluster_statistics(smld, stats_cfg)

# Edge / membrane / interior classification
(smld_passthrough, edge_info) = classify_emitters(smld, KdeValleyConfig())
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

`classify_emitters()` is likewise pass-through: the per-emitter class is mirrored
into `smld.metadata["edge_classify_class"]`, and `info::EdgeClassifyInfo` carries
`class::Vector{Symbol}` (`:outside` / `:membrane` / `:interior`) plus the boundary
geometry. The concrete config type selects the strategy by dispatch. See the
[edge-classification docs](docs/src/edge_classify.md).

```julia
# Cluster boundary polygons (alpha-shape), from a cluster()-labeled SMLD
(smld_out, binfo) = boundary_clusters(smld_out, BoundaryClustersConfig())

# Per-cluster shape & spacing statistics derived from the boundary
(smld_out, sinfo) = boundary_statistics(smld_out, binfo)

# Visualize: Gaussian render + optional centroid/boundary overlay
fig = plot_clusters(smld, smld_out, ClusterPlotConfig(boundaries = binfo))
```

`boundary_clusters()` is pass-through like `cluster_statistics()`: it returns
the same `smld_out` reference alongside a `ClusterBoundaryInfo`.
`boundary_statistics()` is a plain function — not config-dispatched — that
consumes that `ClusterBoundaryInfo` and returns a
`ClusterBoundaryStatisticsInfo`. `plot_clusters()` renders `smld` and overlays
centroids/boundaries sourced from `smld_out`/`binfo`.

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

### MRF density-regime

Adaptive-density clustering for data with multiple density regimes (e.g.
tight ~25 nm aggregates coexisting with μm-scale extended structure).
Avoids the single-ε problem by inferring per-emitter regime labels from
local density via a Gaussian-mixture fit, then smoothing labels over the
Delaunay (or k-NN) neighbor graph using a multi-class Potts MRF before
extracting clusters via connected components on the foreground.

```julia
# Default 2-regime auto-tuning: GMM finds the foreground/background split
# per dataset; smoothness λ is auto-set; CC over Delaunay neighbors.
cfg = MRFDensityClusterConfig(min_points = 10)
(smld_out, info) = cluster(smld, cfg)
regimes = smld_out.metadata["mrf_regime_per_emitter"]   # 0..n_regimes per emitter

# 3-regime with manual thresholds (e.g. learned from training data).
cfg = MRFDensityClusterConfig(
    n_regimes         = 3,
    regime_thresholds = [3.5, 5.0],   # log-density splits, length n_regimes - 1
    min_points        = 10,
)

# Calibrated soft emissions: fit on a calibration ROI, then apply to queries.
# Unlike hard thresholds, borderline interior points can still be rescued by
# high-density neighbors through the MRF smoothness term.
gaussians = calibrate_regime_gaussians(calibration_smld;
                                       n_regimes = 2,
                                       density_estimator = :knn,
                                       density_k = 20)
cfg = MRFDensityClusterConfig(
    density_estimator = :knn,
    density_k         = 20,
    regime_gaussians  = gaussians,
)
```

Lowest regime is treated as background/noise; foreground = regime ≥ 2 is
fed to connected components. Per-dataset by default — each cell gets its
own GMM fit, so the algorithm auto-adapts to whatever density scale that
cell happens to live at.

**2D only.** `use_3d = true` raises `ArgumentError`.
Groups with fewer than 3 points are tagged all-noise.
Groups containing exact-duplicate (x,y) coordinates raise `ArgumentError`.

Output metadata: `mrf_regime_per_emitter` (per-emitter regime ID, original
emitter order), `mrf_lambda_used` (per-group smoothness weight),
`mrf_regime_means` (per-group Gaussian component means in log-density space).

#### When it works

Characterized on a controlled 5×5 μm A431-mimic synthetic (kNN density
estimator, k=20, 2 regimes, ~13–22 k emitters across 12 patches, AR 1–20).
Headline accuracy by density ratio (high / low):

| ratio | kNN-MRF | voronoi-GMM (no MRF) |
|-------|---------|----------------------|
| 1.2×  | 35%     | 65%                  |
| 1.5×  | 69%     | 67%                  |
| 2.0×  | 89%     | 72%                  |
| 3.0×  | 95%     | 69%                  |
| 5.0×  | 96%     | 87%                  |

**Use kNN-MRF when the density ratio is ≥ 2×.** Operational floor: ratio
≥ 1.65× clears 75% accuracy; ratio ≥ 1.85× clears 85%.

**Use calibrated soft emissions (`calibrate_regime_gaussians` +
`regime_gaussians`) when calibration ROIs are available and the density ratio
is below 2×.** This skips per-ROI GMM degeneracy while keeping soft unary
costs, so the MRF can still fill borderline interior points. If no calibration
ROI exists, use voronoi-GMM (`VoronoiDensityConfig` + your own GMM split, or
external thresholding on the per-emitter density extra) because it degrades
more gracefully than a low-contrast auto-MRF collapse.

**Use the kNN density estimator (`density_estimator = :knn, density_k = 20`)
for elongated patches with widths comparable to the local nearest-neighbor
distance** (thin fibers, narrow channels). Round 012 result on the
2× synthetic: kNN closes a 22.5% → 6.5% patch-interior FN-rate gap that
Voronoi-density leaves open. Bound: kNN ball radius (≈ √(k / π · ρ))
must fit inside the structure half-width. For very thin patches drop to
k = 8–10. Voronoi remains the default for blob-shaped clusters where
boundary spillage is a non-issue.

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

## Boundary extraction & shape statistics

Alpha-shape boundary polygons for each cluster, plus per-cluster geometry and
spacing statistics. Built on top of `cluster()`'s output — not a `cluster()`/
`cluster_statistics()` backend itself.

### `boundary_cluster` (low-level)

Analogous to MATLAB's `boundary(x, y, shrink)`. Builds a Delaunay
triangulation and keeps triangles below a `shrink`-controlled circumradius
threshold, then traces the resulting boundary edges into a closed polygon.

```julia
k = boundary_cluster(x, y; shrink = 0.5)   # shrink=0 → convex hull, 1 → tightest
# x[k], y[k] traces the boundary polygon (closed: k[1] == k[end])
```

### `boundary_clusters` / `ClusterBoundaryInfo`

Applies `boundary_cluster` to every cluster in a labeled SMLD.

```julia
(smld_out, _) = cluster(smld, DBSCANConfig(eps_nm = 50.0, min_points = 10))
cfg_bnd       = BoundaryClustersConfig(shrink = 0.5, per_dataset = false)
(_, binfo)    = boundary_clusters(smld_out, cfg_bnd)
# binfo.boundaries_x[j] / binfo.boundaries_y[j]: closed polygon (μm) for
# cluster binfo.cluster_keys[j]
```

Pass-through like `cluster_statistics`: `smld_out` is returned unmodified.
`per_dataset` must match the value used in the preceding `cluster()` call.

### `boundary_statistics` / `ClusterBoundaryStatisticsInfo`

```julia
(_, sinfo) = boundary_statistics(smld_out, binfo;
                                  per_dataset = false,
                                  algorithm   = info.algorithm)
sinfo.area[j]         # μm², at binfo.shrink
sinfo.compactness[j]  # 4π·area/perimeter², 1 for a circle
sinfo.center[j]       # (x, y) centroid, μm
```

Key field groups on `ClusterBoundaryStatisticsInfo` (one entry per cluster
unless noted):

| Group | Fields |
|-------|--------|
| Area / perimeter (μm², μm) | `area`, `perimeter` (at `binfo.shrink`); `areaConvex`/`perimeterConvex` (shrink=0, convex hull); `areaTight`/`perimeterTight` (shrink=1, tightest alpha shape) |
| Shape descriptors | `equiv_radius`, `compactness`, `circularity`, `convexity`, `solidity` |
| Point counts & spacing | `center`, `n_pts_per_cluster`, `n_pts_per_area`, `sigma_actual`, `cluster_width`, `nn_within_clusters` |
| Inter-cluster distances (μm) | `min_c2c_dists` / `min_e2e_dists` (per cluster); `min_c2c_dist` / `min_e2e_dist` (overall scalars) |
| Dataset summary | `n_locs_in`, `n_clustered`, `n_noise`, `algorithm`, `shrink` — duplicated from `ClusterInfo` so the result is a self-contained report |

Fields that need a cluster's raw member points are `NaN`/empty for clusters
with no matching emitters in `smld_out` rather than erroring;
`area`/`perimeter`/`equiv_radius`/`compactness` need only `binfo`'s polygon
and are always computed. See the `ClusterBoundaryStatisticsInfo` docstring
for the full field-by-field reference.

## Visualization

`plot_clusters(smld, smld_out, cfg::ClusterPlotConfig) -> Figure` renders a
Gaussian super-resolution image of `smld` and overlays cluster centroids
and/or boundary polygons from `smld_out` / a `ClusterBoundaryInfo`.

```julia
cfg = ClusterPlotConfig(
    pixel_size_nm  = 10.0,
    color_by       = :id,     # color localizations by cluster id
    categorical    = true,
    boundaries     = binfo,   # optional: draw each cluster's boundary polygon
    show_centroids = true,
    filename       = "clusters.png",
)
fig = plot_clusters(smld, smld_out, cfg)
```

Uses `SMLMRender` for the background render and `CairoMakie` for the figure
and overlays; if `cfg.filename` is set, the figure is also saved to disk.

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
- [DelaunayTriangulation.jl](https://github.com/JuliaGeometry/DelaunayTriangulation.jl) — Voronoi tessellation (Voronoi/MRF backends)
- [AdaptivePredicates.jl](https://github.com/JuliaGeometry/AdaptivePredicates.jl) — exact geometric predicates for EdgeClassify's built-in alpha-shape Delaunay
- [NearestNeighbors.jl](https://github.com/KristofferC/NearestNeighbors.jl) — KDTree NN queries (Hopkins backend, boundary-statistics nearest-neighbor spacing)
- [CairoMakie.jl](https://github.com/MakieOrg/Makie.jl) — figure rendering for `plot_clusters`
- [SMLMRender.jl](https://github.com/JuliaSMLM/SMLMRender.jl) — Gaussian super-resolution background render for `plot_clusters`
