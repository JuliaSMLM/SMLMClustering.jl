# SMLMClustering — API Overview

LLM-parseable API reference. Kept in sync with source by hand; authoritative
source is the docstrings in `src/`. See `README.md` for user-facing narrative.

---

## Entry point

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
| `algorithm` | `Symbol` | `:dbscan`, `:voronoi`, or `:hierarchical` |
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
| `DelaunayTriangulation` | Voronoi tessellation and Delaunay adjacency |

---

## Not yet implemented

- **HDBSCANConfig** — slot reserved; no lightweight pure-Julia HDBSCAN library exists
  as of 2026-04-17. See `KNOWLEDGE_BASE.md` D1 for details. Monitor
  `Clustering.jl` issue #139.
