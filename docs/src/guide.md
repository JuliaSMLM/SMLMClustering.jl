```@meta
CurrentModule = SMLMClustering
```

# User Guide

This guide covers the calling conventions shared by every backend: the inputs, the
three verbs and what they return, the configuration fields common to all labeling
backends, and how to sanity-check a result. For the algorithms themselves see the
[Methods overview](@ref "Methods overview").

## Inputs and units

Every verb takes an `SMLMData.BasicSMLD` as its first argument. Emitter coordinates
are in **µm**. Parameters that express a physical length are given in **nm**
(e.g. `eps_nm`, `cut_threshold` for distance linkages) and converted internally;
each backend page states the unit of every field.

An SMLD may hold multiple **datasets** (cells, ROIs, acquisitions). By default each
backend processes datasets independently — see [`per_dataset`](@ref
"Shared configuration fields") below.

## The three verbs

The three verbs share one calling convention — `verb(smld, cfg) → (smld, info)`,
dispatched on the concrete config type — but do different things with the result.
`cluster` writes an **integer instance label** onto `emitter.id`; `classify_emitters`
leaves `emitter.id` alone and stores a **fixed semantic class** in `metadata`;
`cluster_statistics` writes nothing. Because labeling and classification use different
fields, they **compose** (classify → filter → cluster, or cluster → inspect by
region) — see *Labeling vs. classification* in the
[Methods overview](@ref "Methods overview").

### `cluster` — labeling

```julia
smld_out, info = cluster(smld, cfg)
```

`cluster` is **non-mutating**: the input emitters are deep-copied, cluster labels are
written onto the copy's `emitter.id` (`0` = noise, `1..K` = clusters), and
`info::ClusterInfo` carries the summary. When `cfg.remove_unclustered = true` the
returned `smld_out` contains only clustered emitters.

`ClusterInfo` fields:

| Field | Type | Meaning |
|-------|------|---------|
| `n_locs_in` | `Int` | input localization count |
| `n_clustered` | `Int` | localizations assigned to a cluster (`id > 0`) |
| `n_noise` | `Int` | noise localizations (`id == 0`) |
| `n_clusters` | `Int` | number of distinct clusters |
| `cluster_sizes` | `Vector{Int}` | size of each cluster, indexed by cluster id |
| `algorithm` | `Symbol` | `:dbscan`, `:voronoi`, `:hierarchical`, … |
| `elapsed_s` | `Float64` | wall-clock time of the `cluster` call (s) |

### `cluster_statistics` — read-only statistics

```julia
smld, info = cluster_statistics(smld, stats_cfg)
```

**Pass-through**: the first return value is the same SMLD reference as the input (no
allocation, no mutation); the two-tuple shape is kept for ecosystem symmetry.
`info::ClusterStatisticsInfo` carries:

| Field | Type | Meaning |
|-------|------|---------|
| `n_locs_in` | `Int` | input localization count |
| `statistic` | `Float64` | primary scalar result (e.g. Hopkins `H`) |
| `statistic_name` | `Symbol` | identifier for `statistic` |
| `algorithm` | `Symbol` | backend identifier |
| `elapsed_s` | `Float64` | wall-clock time (s) |
| `extras` | `Dict{Symbol,Any}` | per-backend supplementary outputs (vectors, per-group breakdowns) |

**Convention for vector-valued backends:** a meaningful summary scalar (mean,
median, …) goes in `statistic`; the full per-emitter / per-group vector goes in
`extras` under a descriptive key. This keeps `info.statistic` ergonomic while
preserving the full result.

### `classify_emitters` — edge / membrane / interior

```julia
smld, info = classify_emitters(smld, cfg)   # cfg :: AbstractEdgeClassifyConfig
```

Pass-through; the per-emitter class is mirrored into
`smld.metadata["edge_classify_class"]` and `info::EdgeClassifyInfo` carries
`class::Vector{Symbol}` plus the boundary geometry. See
[Edge / Membrane Classification](@ref).

## Shared configuration fields

Every labeling backend config carries these fields with the same defaults:

| Field | Default | Meaning |
|-------|---------|---------|
| `min_points` | `5` | minimum points for a valid cluster |
| `use_3d` | `false` | include the z-coordinate |
| `per_dataset` | `true` | cluster within each dataset independently |
| `remove_unclustered` | `false` | drop noise emitters from the returned SMLD |

When `per_dataset = true`, `(dataset, id)` uniquely identifies a cluster across a
multi-dataset SMLD; ids are local to each dataset.

!!! note "Backend-specific defaults"
    A few backends override these defaults where the algorithm warrants it — e.g.
    [Point hysteresis](@ref "Point hysteresis") defaults to `min_points = 100` and
    `per_dataset = false`, and [HDBSCAN](@ref) reinterprets `min_points` as the
    core-distance *k*. Each method page lists its own defaults.

## Is the result sane?

A quick checklist after a `cluster` run:

- **Noise fraction.** `info.n_noise / info.n_locs_in` near 1.0 means the length
  scale is too tight (or the data really is unclustered — confirm with the
  [Hopkins statistic](@ref)); near 0.0 with few, huge clusters means it is too loose.
- **Cluster-size distribution.** `info.cluster_sizes` should not be dominated by a
  single giant cluster that swallowed the field (a classic single-`ε` failure — see
  [MRF density-regime](@ref "MRF density-regime")).
- **Count vs. expectation.** `info.n_clusters` should be the right order of magnitude
  for the structure you expect.
- **3D.** Several backends are **2D only** and raise `ArgumentError` on
  `use_3d = true`; check the method page before enabling it.

## Large datasets

- [DBSCAN](@ref), [Voronoi (SR-Tesseler)](@ref "Voronoi (SR-Tesseler)") and
  [Point hysteresis](@ref "Point hysteresis") avoid the O(*n*²) distance matrix and
  are the right choice for ≫10,000 localizations per group.
- [Hierarchical](@ref) builds a dense pairwise distance matrix **per group** —
  prefer it only for small groups.
- With `per_dataset = true`, per-group cost is what matters, not the global *n*.
