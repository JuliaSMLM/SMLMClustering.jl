# SMLMClustering — API Overview

LLM-parseable API reference. Kept in sync with source by hand; authoritative
source is the docstrings in `src/`. See `README.md` for user-facing narrative.

---

## Entry points

The package exposes two parallel entry points:
- **`cluster(smld, cfg)`** — labeling backends (DBSCAN, Hierarchical, Voronoi). Writes per-emitter cluster labels onto a deep-copied SMLD.
- **`cluster_statistics(smld, cfg)`** — read-only spatial-statistic backends (Hopkins, ...). Returns the input SMLD reference unchanged.

See "Sibling entry point — `cluster_statistics`" below for the full asymmetry table.

### `cluster(smld, cfg) -> (smld_out, info)`

**Signature:** `cluster(smld::SMLMData.BasicSMLD, cfg::AbstractClusterConfig) -> (SMLMData.BasicSMLD, ClusterInfo)`

**What it does:** Clusters the localizations in `smld` using the algorithm selected
by the concrete type of `cfg`. Writes per-emitter cluster labels to `emitter.id`
on the **output** SMLD (`0` = noise, `1..K` = cluster, local to each dataset when
`per_dataset=true`). Returns a tuple of the output SMLD and a `ClusterInfo` summary.

**Side effects:** None. The input `smld` is not modified — each backend deep-copies
the input emitters at entry and writes labels onto the copy. When
`cfg.remove_unclustered = true` the returned SMLD contains only the clustered
emitters from that copy; otherwise it contains all emitters from the copy.

**Dispatch:** Concrete config types dispatch to their backend. Passing an unsupported
`AbstractClusterConfig` subtype raises an error naming the supported backends.

---

## Config types

### `AbstractClusterConfig <: SMLMData.AbstractSMLMConfig`

Abstract supertype. Every backend defines a concrete `*Config` subtype.

Shared fields expected on every concrete subtype:

| Field | Type | Default | Meaning |
|-------|------|---------|---------|
| `min_points` | `Int` | `5` | Minimum cluster size (noise threshold) |
| `use_3d` | `Bool` | `false` | Include z-coordinate |
| `per_dataset` | `Bool` | `true` | Cluster within each dataset independently |
| `remove_unclustered` | `Bool` | `false` | Drop noise emitters from output |

---

### `DBSCANConfig <: AbstractClusterConfig`

**Source:** `src/backends/dbscan.jl`
**Algorithm:** DBSCAN via `Clustering.dbscan`
**Library:** Clustering.jl

**Fields:**

| Field | Type | Default | Meaning |
|-------|------|---------|---------|
| `eps_nm` | `Float64` | (required) | Neighborhood radius in **nm**; converted to μm internally |
| `min_points` | `Int` | `5` | Core-point threshold and minimum cluster size |
| `use_3d` | `Bool` | `false` | 2D (`x,y`) or 3D (`x,y,z`) clustering |
| `per_dataset` | `Bool` | `true` | Per-dataset namespacing |
| `remove_unclustered` | `Bool` | `false` | Drop noise emitters |

**Validation:** `eps_nm > 0`, `min_points >= 1`.

**Scalability:** O(n log n) with a KD-tree; suitable for large datasets.

**Supports 3D:** yes.

**Constructor:**
```julia
DBSCANConfig(eps_nm=50.0)
DBSCANConfig(eps_nm=100.0, min_points=3, use_3d=true)
```

---

### `HDBSCANConfig <: AbstractClusterConfig`

**Source:** `src/backends/hdbscan.jl`
**Algorithm:** HDBSCAN* (hierarchical density clustering): core / mutual-reachability distances → MST → condensed tree → stability-based flat extraction
**Library:** pure-Julia implementation (NearestNeighbors.jl for the kNN)

**Fields:**

| Field | Type | Default | Meaning |
|-------|------|---------|---------|
| `min_points` | `Int` | `5` | core-distance neighbor count *k* (density-smoothing scale) |
| `min_cluster_size` | `Int` or `Nothing` | `nothing` | minimum cluster size in the condensed tree; falls back to `min_points` when `nothing` |
| `knn_graph_k` | `Int` | `30` | neighbors per point in the sparse kNN graph used to build the MST |
| `cluster_selection_method` | `Symbol` | `:eom` | flat extraction: `:eom` (excess of mass) or `:leaf` |
| `allow_single_cluster` | `Bool` | `false` | allow the root to be selected as a single cluster |
| `halo_trim_frac` | `Float64` | `0.10` | trim weakly-attached halo points (fell near a cluster's birth) to noise so members track the cluster's physical extent; `0` disables (raw HDBSCAN* subtree labels) |
| `use_3d` | `Bool` | `false` | 2D (`x,y`) or 3D (`x,y,z`) clustering |
| `per_dataset` | `Bool` | `true` | Per-dataset namespacing |
| `remove_unclustered` | `Bool` | `false` | Drop noise emitters |

**Validation:** `min_points >= 1`, `knn_graph_k >= 1`, `cluster_selection_method ∈ {:eom, :leaf}`, and an effective `min_cluster_size >= 2` (after the `nothing → min_points` fallback). A too-small kNN graph self-repairs by bridging disconnected components.

**Returned metadata:** `smld_out.metadata["hdbscan_cluster_persistence"]` and `["hdbscan_cluster_lambda_birth"]` — per-cluster, concatenated across datasets in `per_dataset` order.

**Supports 3D:** yes.

**Constructor:**
```julia
HDBSCANConfig()
HDBSCANConfig(min_points=8, min_cluster_size=20)
```

---

### `HierarchicalConfig <: AbstractClusterConfig`

**Source:** `src/backends/hierarchical.jl`
**Algorithm:** Agglomerative hierarchical clustering via `Clustering.hclust` + `Clustering.cutree`
**Library:** Clustering.jl, Distances.jl

**Fields:**

| Field | Type | Default | Meaning |
|-------|------|---------|---------|
| `cut_threshold` | `Union{Float64,Nothing}` | `nothing` | Dendrogram cut height; **unit depends on linkage** (see caveat) |
| `n_clusters` | `Union{Int,Nothing}` | `nothing` | Cut into exactly K clusters; mutually exclusive with `cut_threshold` |
| `linkage` | `Symbol` | `:ward` | `:single`, `:complete`, `:average`, or `:ward` |
| `min_points` | `Int` | `5` | Sub-threshold clusters relabeled noise |
| `use_3d` | `Bool` | `false` | 2D or 3D clustering |
| `per_dataset` | `Bool` | `true` | Per-dataset namespacing |
| `remove_unclustered` | `Bool` | `false` | Drop noise emitters |

**Validation:** exactly one of `cut_threshold` / `n_clusters` set (both-or-neither →
`ArgumentError`), `cut_threshold > 0` when set, `n_clusters >= 1` when set,
`min_points >= 1`, `linkage` in `(:single, :complete, :average, :ward)`.

**Scalability:** O(n²) pairwise distance matrix per group. Prefer `DBSCANConfig` for
groups with ≫10,000 localizations.

**Supports 3D:** yes.

**`cut_threshold` unit convention:** for distance-based linkages (`:single`,
`:complete`, `:average`) the value is in **nanometers** and is converted to μm
internally (`h = cut_threshold / 1000.0`). For `:ward` the dendrogram height is a
variance-increase cost (roughly μm²) and is passed through without conversion —
there is no meaningful nm interpretation under Ward, which is why `n_clusters` is
usually the cleaner choice for Ward.

**Constructor:**
```julia
HierarchicalConfig(cut_threshold=200.0, linkage=:single)  # distance-based, nm
HierarchicalConfig(n_clusters=3, linkage=:ward)           # Ward, by count
```

---

### `MRFDensityClusterConfig <: AbstractClusterConfig`

**Source:** `src/backends/mrf_density.jl`
**Algorithm:** local density → soft regime unaries → multi-class Potts MRF (ICM) → CC on foreground
**Library:** DelaunayTriangulation.jl, NearestNeighbors.jl, Statistics

Adaptive-density clustering for data with multiple density regimes (e.g.
tight ~25 nm aggregates next to μm-scale extended structure). The MRF
smoothness term enforces spatial coherence: borderline middle points stay
foreground, isolated tight knots in a low-density sea get demoted to
background. No global ε, no per-dataset density tuning.

**Fields:**

| Field | Type | Default | Meaning |
|-------|------|---------|---------|
| `n_regimes` | `Int` | `2` | Number of density regimes (lowest = noise/background) |
| `regime_thresholds` | `Union{Nothing, Vector{Float64}}` | `nothing` | Optional explicit log-density thresholds (length `n_regimes - 1`, sorted asc); hard-bins labels and bypasses GMM |
| `regime_gaussians` | `Union{Nothing, NamedTuple}` | `nothing` | Optional calibrated `(means, vars, weights)` Gaussian emissions in log-density space; bypasses GMM while preserving soft unaries |
| `density_estimator` | `Symbol` | `:voronoi` | Per-emitter density estimator: `:voronoi` or `:knn` |
| `density_k` | `Int` | `20` | k for `density_estimator=:knn` |
| `smoothness_lambda` | `Union{Nothing, Float64}` | `nothing` | MRF smoothness weight; when `nothing`, auto-tuned per group via MAD of unary range |
| `graph_kind` | `Symbol` | `:delaunay` | Neighbor graph: `:delaunay` (free, reuses tessellation) or `:knn` |
| `graph_k` | `Int` | `8` | k for kNN graph (used only with `graph_kind=:knn`) |
| `inference` | `Symbol` | `:icm` | MRF inference; only `:icm` supported in v1 |
| `icm_iters` | `Int` | `50` | Maximum ICM passes (early termination on convergence) |
| `min_points` | `Int` | `5` | Minimum cluster size after CC |
| `use_3d` | `Bool` | `false` | **Must be `false`** — 3D not supported |
| `per_dataset` | `Bool` | `true` | Run pipeline per dataset (independent GMM fit per cell) |
| `remove_unclustered` | `Bool` | `false` | Drop noise emitters |

**Validation:** `n_regimes >= 2`, at most one of `regime_thresholds` /
`regime_gaussians` set, `regime_thresholds` (when set) length =
`n_regimes - 1` and sorted ascending, `regime_gaussians` (when set) has sorted
finite means plus positive variances/weights, `density_estimator in
(:voronoi, :knn)`, `density_k >= 1`, `graph_kind in (:delaunay, :knn)`,
`inference === :icm`, `graph_k >= 1`, `icm_iters >= 1`, `min_points >= 1`,
`smoothness_lambda > 0` when set, `use_3d == false`.

**Scalability:** dominated by Voronoi tessellation (O(n log n)) and ICM
passes (O(`icm_iters` × n × `n_regimes` × avg-degree)). Comparable to
DBSCAN in practice for n up to 100k per group.

**Supports 3D:** no. `use_3d = true` raises `ArgumentError` directing users
to `DBSCANConfig` or `HierarchicalConfig`.

**Output metadata** (in `smld_out.metadata`):

| Key | Type | Meaning |
|-----|------|---------|
| `"mrf_regime_per_emitter"` | `Vector{Int}` | Per-emitter regime 0..`n_regimes` (0 = ungroupable, 1 = lowest density, `n_regimes` = highest), in original emitter order |
| `"mrf_lambda_used"` | `Vector{Float64}` | Per-group λ actually used (auto or explicit), in `_group_by_dataset` order |
| `"mrf_regime_means"` | `Vector{Vector{Float64}}` | Per-group Gaussian means (sorted ascending, log-density space). When hard `regime_thresholds` were provided, the per-group entry is filled with `NaN`s |

**Constructor:**
```julia
MRFDensityClusterConfig()                                          # 2-regime, GMM auto
MRFDensityClusterConfig(n_regimes = 3, min_points = 10)
MRFDensityClusterConfig(n_regimes = 3, regime_thresholds = [3.5, 5.0])
gaussians = calibrate_regime_gaussians(cal_smld; density_estimator = :knn)
MRFDensityClusterConfig(density_estimator = :knn, regime_gaussians = gaussians)
MRFDensityClusterConfig(graph_kind = :knn, graph_k = 12, smoothness_lambda = 0.5)
```

---

### `VoronoiConfig <: AbstractClusterConfig`

**Source:** `src/backends/voronoi.jl`
**Algorithm:** Voronoi-tessellation density clustering (SR-Tesseler)
**Library:** DelaunayTriangulation.jl
**Reference:** Levet et al., Nat. Methods 2015

**Fields:**

| Field | Type | Default | Meaning |
|-------|------|---------|---------|
| `density_factor` | `Float64` | `2.0` | Density threshold multiplier (dense ⟺ area < mean/factor) |
| `min_points` | `Int` | `5` | Minimum cluster size |
| `use_3d` | `Bool` | `false` | **Must be `false`** — 3D not supported |
| `per_dataset` | `Bool` | `true` | Per-dataset namespacing |
| `remove_unclustered` | `Bool` | `false` | Drop noise emitters |

**Validation:** `density_factor > 0`, `min_points >= 1`, `use_3d == false` (raises
`ArgumentError` otherwise).

**Supports 3D:** no. `use_3d = true` raises `ArgumentError` directing users to
`DBSCANConfig` or `HierarchicalConfig`.

**Degenerate input handling:**
- Groups with fewer than 3 points: all emitters tagged noise (no tessellation).
- Groups with exact-duplicate (x,y) coordinates: `ArgumentError` raised before
  triangulation (deduplicate input first).

**Boundary note:** Cells are clipped to the convex hull; hull generators have
smaller-than-true cell areas, which can slightly bias mean-area estimates on small groups.

**Constructor:**
```julia
VoronoiConfig()
VoronoiConfig(density_factor=3.0, min_points=10)
```

---

### `PointHysteresisConfig <: AbstractClusterConfig`

**Source:** `src/backends/point_hysteresis.jl`
**Algorithm:** seed-and-grow connected components on a directed kNN graph, gated by caller-supplied `seed` / `support` masks (density hysteresis)
**Library:** NearestNeighbors.jl

Unlike the other labeling backends, `cluster` for `PointHysteresisConfig` takes two
**required keyword masks**, computed upstream (e.g. from `LocalContrastFeature`):

```julia
cluster(smld, PointHysteresisConfig(); seed=seed_mask, support=support_mask)
```

A `support`-connected component becomes a cluster iff it contains at least one `seed`
point and has at least `min_points` members.

**Fields:**

| Field | Type | Default | Meaning |
|-------|------|---------|---------|
| `graph_k` | `Int` | `12` | neighbors per point in the directed kNN growth graph |
| `min_points` | `Int` | `100` | minimum component size to keep as a cluster |
| `use_3d` | `Bool` | `false` | 2D or 3D |
| `per_dataset` | `Bool` | `false` | Per-dataset namespacing |
| `remove_unclustered` | `Bool` | `false` | Drop noise emitters |

**Required keyword args (on `cluster`):** `seed::AbstractVector{Bool}` and `support::AbstractVector{Bool}`, each of length = number of emitters.

**Supports 3D:** yes.

**Constructor:**
```julia
PointHysteresisConfig(graph_k=12, min_points=50)
```

---

## Result type

### `ClusterInfo <: SMLMData.AbstractSMLMInfo`

**Source:** `src/types.jl`

**Fields:**

| Field | Type | Meaning |
|-------|------|---------|
| `n_locs_in` | `Int` | Input localization count |
| `n_clustered` | `Int` | Localizations with `id > 0` |
| `n_noise` | `Int` | Localizations with `id == 0` |
| `n_clusters` | `Int` | Distinct clusters formed |
| `cluster_sizes` | `Vector{Int}` | Size of cluster `k` at index `k`; length = `n_clusters` |
| `algorithm` | `Symbol` | `:dbscan`, `:hdbscan`, `:hierarchical`, `:voronoi`, `:mrf_density`, or `:point_hysteresis` |
| `elapsed_s` | `Float64` | Wall-clock time of the `cluster` call (seconds) |

**`cluster_sizes` convention:** When `per_dataset = true`, sizes from all datasets are
concatenated in dataset-visit order. `cluster_sizes[k]` is the size of the k-th cluster
in that order; the mapping to `(dataset, id)` is not stored in `ClusterInfo` — reconstruct
from emitter `dataset` and `id` fields if needed.

**`show`:** `ClusterInfo(n_clustered/n_locs_in clustered, n_clusters clusters, algorithm=:X, T ms)`

---

## Label convention

- `emitter.id == 0` — noise / rejected
- `emitter.id == k` (k ≥ 1) — cluster k, local to the dataset

When `per_dataset = true`, `(emitter.dataset, emitter.id)` uniquely identifies a
cluster across a multi-dataset SMLD. The same `id` in different datasets refers to
different clusters. Cluster step writes labels into `emitter.id` and is designed to
run after FrameConnect / BaGoL (which use `track_id`, leaving `id` free).

---

## Package dependencies

| Package | Role |
|---------|------|
| `SMLMData` | `BasicSMLD`, emitter types, `AbstractSMLMConfig`, `AbstractSMLMInfo` |
| `Clustering` | `dbscan`, `hclust`, `cutree` |
| `Distances` | `pairwise(Euclidean(), X; dims=2)` for hierarchical distance matrix |
| `DelaunayTriangulation` | Voronoi tessellation and Delaunay adjacency (Voronoi/MRF backends) |
| `AdaptivePredicates` | Exact geometric predicates for EdgeClassify's built-in alpha-shape Delaunay |
| `NearestNeighbors` | KDTree NN queries for the Hopkins backend |
| `Random` | Seeded RNG (`Xoshiro`) for reproducible Hopkins repeats |
| `Statistics` | `median` for the Voronoi-density summary statistic |

---

## Sibling entry point — `cluster_statistics`

`cluster_statistics(smld, cfg) -> (smld, ClusterStatisticsInfo)` is a **read-only**
sibling to `cluster()` for spatial-statistic backends (clustering tendency, density
diagnostics, etc.). The two interfaces are intentionally asymmetric:

| Aspect | `cluster()` | `cluster_statistics()` |
|--------|-------------|------------------------|
| Writes labels? | Yes (`emitter.id`) | No |
| Input SMLD | Deep-copied | Pass-through (same reference) |
| Result info | `ClusterInfo` | `ClusterStatisticsInfo` |
| Abstract supertype | `AbstractClusterConfig` | `AbstractStatisticsConfig` |

Both return a `(smld, info)` tuple for ecosystem symmetry, but the SMLD returned
by `cluster_statistics` is `===` the input — no allocation, no mutation.

### `AbstractStatisticsConfig <: SMLMData.AbstractSMLMConfig`

Abstract supertype. Every spatial-statistic backend defines a concrete `*Config` subtype.

Shared fields expected on every concrete subtype:

| Field | Type | Default | Meaning |
|-------|------|---------|---------|
| `use_3d` | `Bool` | `false` | Include z-coordinate |
| `per_dataset` | `Bool` | `true` | Compute per dataset and aggregate |

### `ClusterStatisticsInfo <: SMLMData.AbstractSMLMInfo`

**Source:** `src/types.jl`

**Fields:**

| Field | Type | Meaning |
|-------|------|---------|
| `n_locs_in` | `Int` | Input localization count |
| `statistic` | `Float64` | Primary scalar result |
| `statistic_name` | `Symbol` | Identifier for `statistic` (`:hopkins`, ...) |
| `algorithm` | `Symbol` | Backend identifier (`:hopkins`, ...) |
| `elapsed_s` | `Float64` | Wall-clock time (seconds) |
| `extras` | `Dict{Symbol,Any}` | Vector / supplementary outputs |

**Convention for vector-valued backends:** put a meaningful summary scalar in
`statistic` (mean, median, ...) and the full vector in `extras` under a
descriptive key. Hopkins per-dataset vector goes under `:hopkins_per_dataset`.

**`show`:** `ClusterStatisticsInfo(:name=value, n_locs_in=N, algorithm=:X, T ms)`

### `HopkinsConfig <: AbstractStatisticsConfig`

**Source:** `src/backends/hopkins.jl`
**Algorithm:** Hopkins clustering-tendency statistic (sample-based)
**Library:** NearestNeighbors.jl (KDTree)

**Fields:**

| Field | Type | Default | Meaning |
|-------|------|---------|---------|
| `n_samples` | `Int` | `20` | Reference / sampled point count per repeat |
| `random_repeats` | `Int` | `1` | Independent repeats to average |
| `seed` | `Union{Int,Nothing}` | `nothing` | RNG seed for reproducibility |
| `use_3d` | `Bool` | `false` | Include z-coordinate |
| `per_dataset` | `Bool` | `true` | Compute per dataset and aggregate |
| `region` | `Nothing` / `Symbol` / polygon / `Dict{Int,polygon}` | `nothing` | observation window for the uniform reference points (2D). `nothing` = data bbox; a `Vector{NTuple{2,Float64}}` = rejection-sample references inside it; `:metadata` = use `smld.metadata["edge_outer_polygon"]` (from `classify_emitters`); `Dict(dataset_id=>polygon)` = per dataset. Corrects false "clustered" on non-convex domains; incompatible with `use_3d=true` |

**Validation:** `n_samples >= 1`, `random_repeats >= 1`, `n_samples <= n_points` per group
(violations within a group return NaN for that group rather than erroring).

**Returned info:**
- `statistic`: mean H across datasets when `per_dataset=true`, single H when `false`.
- `extras[:hopkins_per_dataset]`: `Vector{Float64}` of per-dataset H (only populated when `per_dataset=true`).

**Interpretation:**
- `H ≈ 0.5` — uniform / Poisson-consistent
- `H → 1.0` — strong clustering tendency
- `H → 0.0` — anti-clustering / regular spacing

**Edge cases:** empty group / `n_samples > n_points` / degenerate (zero-extent)
bbox → group H is `NaN`. Aggregate `statistic` is the mean of non-NaN per-dataset
values; if all groups are NaN, `statistic` is NaN.

**Constructor:**
```julia
HopkinsConfig()
HopkinsConfig(n_samples=50, random_repeats=10, seed=1)
```

### `VoronoiDensityConfig <: AbstractStatisticsConfig`

**Source:** `src/backends/voronoi_density.jl`
**Algorithm:** Per-emitter Voronoi cell area → local density `ρᵢ = 1/Aᵢ`
**Library:** DelaunayTriangulation.jl, Statistics

**Fields:**

| Field | Type | Default | Meaning |
|-------|------|---------|---------|
| `use_3d` | `Bool` | `false` | **Must be `false`** — 3D Voronoi not supported |
| `per_dataset` | `Bool` | `true` | Tessellate each dataset independently |

**Validation:** `use_3d == false` (raises `ArgumentError` directing to `VoronoiConfig` docstring otherwise).

**Returned info:**
- `statistic`: median density across all emitters that received a valid Voronoi cell (μm⁻²).
- `statistic_name`: `:median_density`.
- `algorithm`: `:voronoi_density`.
- `extras[:density_per_emitter]`: `Vector{Float64}` of length `n_locs_in`, in **original emitter order** (NOT grouped by dataset). Units: μm⁻². Emitters in groups smaller than 3 receive `NaN`.
- `extras[:area_per_emitter]`: `Vector{Float64}` of length `n_locs_in`, same ordering. Units: μm².

**Edge cases:**
- Groups with fewer than 3 points: those emitters receive `NaN` density and area; other groups proceed normally.
- Empty SMLD: empty per-emitter vectors, `statistic = NaN`.
- Group with exact-duplicate `(x, y)` coordinates: `ArgumentError` raised before triangulation (mirrors `VoronoiConfig`'s guard).

**Use case:** intended for downstream callers running their own thresholding on the per-emitter density (Otsu / GMM on `log ρ`, fixed cutoff, etc.) — e.g. cell-structure masking on dense membrane regions before per-cell clustering.

**Constructor:**
```julia
VoronoiDensityConfig()
VoronoiDensityConfig(per_dataset=false)
```

### `LocalContrastFeature <: AbstractStatisticsConfig`

**Source:** `src/backends/local_contrast.jl`
**Algorithm:** per-emitter local-density contrast — fine-scale kNN log-density minus the median log-density over a coarser neighborhood (cancels a baseline density gradient)
**Library:** NearestNeighbors.jl, Statistics

**Fields:**

| Field | Type | Default | Meaning |
|-------|------|---------|---------|
| `density_k` | `Int` | `200` | fine-scale neighbor count for the per-point log-density |
| `background_k` | `Int` | `2000` | coarse-scale neighbor count for the local baseline; must be `> density_k` |
| `use_3d` | `Bool` | `false` | 2D or 3D (the density normalization uses the 2D disk area regardless) |
| `per_dataset` | `Bool` | `false` | compute each dataset independently |

**Validation:** `background_k > density_k` and `density_k >= 1` (raises `ArgumentError` otherwise).

**Returned info:**
- `statistic`: median of the finite per-emitter contrasts.
- `statistic_name`: `:median_local_contrast`.
- `algorithm`: `:local_contrast`.
- `extras[:contrast_per_emitter]`: `Vector{Float64}`, length `n_locs_in`, original emitter order.
- `extras[:log_density_per_emitter]`: `Vector{Float64}`, the fine kNN log-density, same ordering.

**Use case:** a threshold-ready foreground feature on data with non-stationary baseline density — e.g. to build the `seed` / `support` masks for `PointHysteresisConfig`.

**Constructor:**
```julia
LocalContrastFeature()
LocalContrastFeature(density_k=20, background_k=200)
```

---

## Sibling entry point — `classify_emitters`

`classify_emitters(smld, cfg) -> (smld, EdgeClassifyInfo)` is a third **read-only**
sibling (alongside `cluster` / `cluster_statistics`) that classifies each emitter as
`:outside` / `:membrane` / `:interior`. The concrete config type selects the strategy
by dispatch. The coordinate form `classify_emitters(x_um, y_um, cfg; fov_um) ->
EdgeClassifyInfo` is the core; the SMLD form mirrors the class into
`smld.metadata["edge_classify_class"]` (`Vector{String}`) and returns the smld
pass-through. Artifact writing is a separate step:
`write_edge_artifacts(leaf, info, x_um, y_um; condition, cell)`.

### `AbstractEdgeClassifyConfig <: SMLMData.AbstractSMLMConfig`

Abstract supertype; each strategy is a concrete subtype dispatched as a
`classify_emitters` method. Shared geometry fields on every subtype:

| Field | Type | Default | Meaning |
|-------|------|---------|---------|
| `alpha_nm` | `Float64` | 300 (kde: 600) | alpha-shape circumradius |
| `membrane_nm` | `Float64` | 100 | membrane band width (nm) |
| `reflect_radius_nm` | `Float64` | 1500 | FOV-reflection band (nm) |
| `fov_trunc_tol_nm` | `Float64` | 150 | FOV-truncation tolerance (nm) |

### `EdgeClassifyInfo <: SMLMData.AbstractSMLMInfo`

**Source:** `src/edge_classify/info.jl`

| Field | Type | Meaning |
|-------|------|---------|
| `class` | `Vector{Symbol}` | authoritative per-emitter class (`:outside`/`:membrane`/`:interior`) |
| `inside_outer` | `BitVector` | **geometric** containment in the *classification* loop `loops[1]` |
| `dist_to_outer_um` | `Vector{Float64}` | distance to that classification loop; `NaN` if not inside |
| `outer_polygon` | polygon | **published** boundary: un-reflected footprint alpha-shape, **FOV-clipped** (drawn boundary + Hopkins `region=:metadata` window) |
| `loops`, `loop_diagnostics` | — | reflected/FOV-augmented loops (`loops[1]` = labeling boundary) + per-loop diagnostics |
| `config` | concrete config | provenance |
| `n_outside` / `n_membrane` / `n_interior` | `Int` | class counts |

Accessors: `in_cell(info)` (= `class .!= :outside`), `interior_fraction(info)`.
**Filter on `class`, never `inside_outer`.** For `KdeValleyConfig` the enclosure stage
folds enclosed background into `:interior` while `inside_outer` stays geometric, so the
enclosure-recovered set is exactly `class == :interior && inside_outer == false`.
**Two boundaries:** `outer_polygon` is the published footprint (un-reflected,
FOV-clipped) used for drawing + the Hopkins window and follows real concavities;
`loops[1]` is the reflected classification loop that `inside_outer`/`dist`/`membrane`
are measured against. They coincide when no FOV side is truncated.

### `OuterPolygonConfig <: AbstractEdgeClassifyConfig`

**Source:** `src/edge_classify/classify.jl`
**Algorithm:** FOV-reflect → multi-K k-NN density gate → alpha-shape outer loop → point-in-polygon + membrane band
**Library:** NearestNeighbors.jl, AdaptivePredicates.jl (built-in pure-Julia Delaunay)

**Extra fields:** `k_list::Tuple{Vararg{Int}}` = `(16, 128)`; `rho_k_thresh::Float64` = `200` (µm⁻²).

**Constructor:**
```julia
OuterPolygonConfig()
OuterPolygonConfig(alpha_nm=400.0, rho_k_thresh=50.0)
```

### `KdeValleyConfig <: AbstractEdgeClassifyConfig`

**Source:** `src/edge_classify/classify.jl` + `src/edge_classify/gates.jl`
**Algorithm:** Gaussian-KDE density → background/cell valley threshold → footprint fill → outer-polygon geometry on the footprint subset → 8-ray enclosure reclass. Validated adaptive dSTORM gate; per-FOV adaptive (no per-cell density tuning).
**Library:** NearestNeighbors.jl, AdaptivePredicates.jl (built-in pure-Julia Delaunay), Statistics

**Extra fields:** `sigma_nm`=150, `rmax_sigma`=3.0, `valley_nbins`=140, `valley_floorfrac`=0.05, `valley_smooth`=4, `footprint_bin_um`=0.2, `footprint_closing_px`=3, `enclosure_bin_um`=0.2, `enclosure_min_hits`=6. Defaults reproduce the validated A431 dSTORM set (note `alpha_nm`=600).

**Constructor:**
```julia
KdeValleyConfig()                       # validated A431 defaults
KdeValleyConfig(sigma_nm=120.0)
```
