# Precision DBSCAN

Precision-weighted DBSCAN is a σ-aware variant of [DBSCAN](@ref): instead of one
global radius, the neighbor test scales with each localization's own uncertainty.
Two localizations are neighbors when

```math
\lVert p_i - p_j \rVert < \texttt{nsigma}\,\bigl(\sigma_{\mathrm{eff},i} + \sigma_{\mathrm{eff},j}\bigr),
```

so precise localizations must be closer to link than imprecise ones. It is a
`cluster` labeling backend selected by a `PrecisionDBSCANConfig`, and — unlike the
Euclidean [DBSCAN](@ref) backend — is a self-contained implementation (the metric is
not Euclidean, so it cannot use `Clustering.dbscan`).

## Concept

DBSCAN's fixed `eps_nm` treats every localization as equally certain. In SMLM the
localization precision `σ` varies point to point (photon count, background), so a
fixed radius is either too generous for precise points or too strict for imprecise
ones. Precision DBSCAN replaces the radius with a per-pair threshold
`nsigma·(σ_eff_i + σ_eff_j)`: the neighborhood breathes with the data's own
uncertainty and there is no absolute length to tune — only the dimensionless
`nsigma`.

Each emitter's scalar `σ_eff` is the geometric mean of its per-axis localization
precisions — `√(σ_x·σ_y)` in 2D, `∛(σ_x·σ_y·σ_z)` in 3D.

## How it works

For a group of localizations the backend:

1. Builds a **neighbor cache** — a KD-tree range query out to
   `max_radius = nsigma · 2 · maxᵢ σ_eff,ᵢ` (a superset that provably contains every
   pair that could pass the threshold), storing each candidate pair `(i, j)` with its
   raw Euclidean distance. The range queries are threaded.
2. **Re-thresholds** the cached pairs with `d < nsigma·(σ_eff_i + σ_eff_j)` to get
   the active edge set.
3. **Labels** the active graph:
   - the **core-point** rule (`min_points ≥ 1`, the config path): a point is a *core
     point* when its active degree is `≥ min_points`; core points sharing an active
     edge merge; a non-core *border* point joins the lowest-id adjacent core cluster;
     a point with no active core neighbor is **noise**;
   - clusters below `min_points` in final size are dropped to noise and the survivors
     are compact-relabeled `1..K`, matching [DBSCAN](@ref).

The label pass is deterministic and independent of thread scheduling.

### The reuse primitive (build once, relabel many)

Steps 1 and 2 are separable, which is the point of the low-level API. The geometry
cache depends only on the coordinates and `max_radius`; `σ_eff` and `nsigma` enter
only at label time. So a caller that repeatedly relabels the **same** points with
different `σ_eff` / `nsigma` (e.g. an EM-style estimator sweeping a scale parameter)
can build the cache **once** and reuse it:

```julia
using SMLMClustering

coords = ...                      # 2×N (or 3×N) matrix, columns = points, same unit as σ_eff
g = build_precision_neighbor_graph(coords, max_radius)   # once; threaded prepass

for nsigma in schedule
    σ_eff = current_precisions()  # may change every iteration
    labels = precision_dbscan_labels(g, σ_eff, nsigma; min_pts = 0)
    # ... use labels ...
end
```

Build `max_radius` from the **coarsest** threshold the loop will ever reach
(`≥ nsigma·2·max(σ_eff)` over all iterations), so the cache stays a valid superset;
`precision_dbscan_labels` asserts this each call (`check_superset = true`) and errors
if the cache is too tight. With `min_pts = 0` the label pass is pure connected
components (union-find) — order-free and bit-identical whether the graph was freshly
built or reused. `precision_dbscan_labels!` writes into a preallocated label vector
for the hot loop.

## Configuration

`PrecisionDBSCANConfig <: AbstractClusterConfig`. Construct with keywords; `nsigma`
is required.

| field | default | unit | meaning |
|---|---|---|---|
| `nsigma` | (required) | — | neighbor radius in units of the summed precision `σ_eff_i + σ_eff_j`. Must be `> 0`. |
| `min_points` | `3` | count | core-point threshold and minimum cluster size (clusters smaller than this become noise), as in [DBSCAN](@ref). Must be `≥ 1`. |
| `use_3d` | `false` | — | cluster in `(x, y, z)` using `σ_z` when `true` (requires 3D emitters), otherwise `(x, y)`. |
| `per_dataset` | `true` | — | cluster within each `dataset` index independently so `(dataset, id)` identifies a cluster. |
| `remove_unclustered` | `false` | — | drop noise emitters (`id == 0`) from the returned SMLD. |

```julia
using SMLMClustering

cfg = PrecisionDBSCANConfig(nsigma = 5.0, min_points = 5)
smld_out, info = cluster(smld, cfg)
println(info)   # ClusterInfo, algorithm = :precision_dbscan
```

## Output & interpretation

`cluster` returns `(smld_out, info)` with the same contract as [DBSCAN](@ref):
`smld_out` is a deep copy carrying per-emitter labels on `emitter.id` (`0` = noise,
`1..K` = cluster; local to each dataset when `per_dataset = true`), and
`info::ClusterInfo` has `algorithm = :precision_dbscan`, `cluster_sizes` in cluster-id
order, and the `n_locs_in` / `n_clustered` / `n_noise` / `n_clusters` counts.

## Notes & caveats

- **Units.** `σ_eff` is derived from the emitters' `σ_x` / `σ_y` (/ `σ_z`), which are
  in **microns** (as are the coordinates); `nsigma` is dimensionless. The low-level
  `build_precision_neighbor_graph` / `precision_dbscan_labels` take `coords`,
  `max_radius`, and `σ_eff` in **one consistent length unit** of your choosing.
- **`min_points` semantics.** As in [DBSCAN](@ref), the same value is the core-point
  threshold and the minimum cluster size. The config path requires `min_points ≥ 1`;
  the primitive additionally supports `min_pts = 0` for pure connected components (no
  noise), which the config wrapper does not expose.
- **Zero precision.** If every localization has `σ_eff = 0`, no neighborhood can form;
  the config raises `ArgumentError`.
- **Parallelism.** The neighbor prepass is threaded (start Julia with `-t` / set
  `JULIA_NUM_THREADS`); the label pass is serial but deterministic.

## References

- M. Ester, H.-P. Kriegel, J. Sander, and X. Xu, "A Density-Based Algorithm for
  Discovering Clusters in Large Spatial Databases with Noise," *KDD-96*, AAAI Press,
  1996, pp. 226–231 — the DBSCAN algorithm this backend generalizes.
