# Round 002 — DBSCAN backend

**Date:** 2026-04-17
**Status:** done
**Priority worked on:** Priority 2 — Implement `DBSCANConfig` backend

## Hypothesis

With the interface from Round 001 in place (`AbstractClusterConfig`, `ClusterInfo`, `cluster(smld, cfg)` generic), adding a concrete DBSCAN backend requires only:

1. A `DBSCANConfig <: AbstractClusterConfig` struct with shared fields (V3) plus algorithm-specific `eps_nm`.
2. One specialized `cluster(smld, ::DBSCANConfig)` method that (a) builds a coordinate matrix in the right dimensionality, (b) runs Clustering.jl's DBSCAN per-dataset (or flat), (c) writes labels to `emitter.id` (V1), (d) populates `ClusterInfo`.

If that is true, the interface is confirmed correct and the three remaining backends (HDBSCAN, Voronoi, hierarchical) are mechanical additions in the same shape.

## What was attempted

- **Project.toml** — added `Clustering = "aaaa29a8-..."` to `[deps]` with `Clustering = "0.15"` compat. Also promoted test-only `Random` into `[extras]`/`[targets]` so `Pkg.test()` sees it (the synthetic-blob tests seed an `Xoshiro` for reproducibility).
- **src/backends/dbscan.jl** (new) — `Base.@kwdef struct DBSCANConfig` with fields `eps_nm::Float64`, `min_points::Int=5`, `use_3d::Bool=false`, `per_dataset::Bool=true`, `remove_unclustered::Bool=false`. Private `_coords_matrix(emitters, use_3d)` helper builds a `d×n` column-major matrix (2×n for 2D, 3×n for 3D) in microns, and raises a clear error if `use_3d=true` but the emitters lack `:z`. The `cluster` method:
  - Validates `eps_nm > 0` and `min_points >= 1` (`ArgumentError`).
  - Converts `eps_nm` (nm) → `radius_μm` (μm) for the KDTree call.
  - Groups indices by `emitter.dataset` when `per_dataset=true` (sorted-keys iteration for deterministic ordering); otherwise treats everything as one group.
  - For each group, runs `Clustering.dbscan(X, radius_μm; min_neighbors=cfg.min_points, min_cluster_size=cfg.min_points)`. Using `min_points` for **both** Clustering.jl parameters matches classical DBSCAN semantics (minPts is both the core-point threshold and the minimum cluster size).
  - Writes `res.assignments[j]` to `emitters[idxs[j]].id` — label namespace stays local per dataset (aligns with V3's "(dataset, id) uniquely identifies" invariant).
  - Appends `res.counts` into the running `cluster_sizes::Vector{Int}`; `length(cluster_sizes) == n_clusters` and `sum(cluster_sizes) == n_clustered` always hold.
  - Rebuilds the output `BasicSMLD` with the same camera/frames/datasets/metadata. If `remove_unclustered=true`, filters to `id != 0` emitters.
- **src/SMLMClustering.jl** — `include("backends/dbscan.jl")` and extended `export` to include `DBSCANConfig`.
- **test/test_dbscan.jl** (new, 45 tests) — testsets:
  1. Config construction: default values and explicit-keyword overrides.
  2. Three well-separated Gaussian blobs (centers 1 μm apart, σ=10 nm, eps=100 nm, min_points=5) + 30 distant noise points → must yield `n_clusters == 3`, `n_clustered >= 115`, `n_noise >= 20`. Checks `algorithm === :dbscan`, `elapsed_s >= 0`, `sum(cluster_sizes) == n_clustered`, `n_clustered + n_noise == n_locs_in`.
  3. `emitter.id` is written in place; `cluster_sizes[k] == count(e -> e.id == k, emitters)`. `remove_unclustered=true` drops noise from the output SMLD and leaves clustered emitters intact.
  4. `per_dataset=true` produces ids `[1, 2]` within each dataset (two datasets with identical spatial layouts produce 4 clusters total, 2 per dataset). Contrast: `per_dataset=false` merges the spatially-identical points across datasets into 2 clusters.
  5. Argument validation: `eps_nm == 0.0` and `min_points == 0` each raise `ArgumentError`.
  6. `use_3d=true` on `Emitter2DFit` data raises `ErrorException` (from `_coords_matrix`).
  7. 3D clustering: two blobs at the same (x, y) but z separated by 1 μm cluster separately under `use_3d=true`.
  8. Empty SMLD: returns `n_clusters=0`, `n_clustered=0`, empty `cluster_sizes`, empty output emitters.
- **test/runtests.jl** — `include("test_dbscan.jl")` inside the top-level `@testset`.

## What worked

- Full suite green: **56/56 passing in 2.2 s via `Pkg.test()`**. The 11 interface tests from Round 001 still pass; the 45 new tests all pass.
- The interface hypothesis is confirmed: the DBSCAN implementation needed no changes to `AbstractClusterConfig`, `ClusterInfo`, or the `cluster` generic from Round 001.
- Label namespacing under `per_dataset=true` works as V3 specifies: two datasets with identical-coordinate blobs produce four clusters with non-overlapping `(dataset, id)` pairs but overlapping `id` values (both datasets use `id ∈ {1, 2}`).
- `remove_unclustered=true` preserves the in-place mutation of input emitters' ids (expected, since the output SMLD shares the emitter vector when the flag is off) but filters the output to `id != 0`.

## What failed

- First test run errored because `sort(...)` was called on a `Base.Generator` inside the per-dataset testset (Julia 1.12 requires a collection, not a generator). Fixed by switching to `sort!(unique(e.id for e in ... if ...))` — `unique` on a generator returns a `Vector`, which `sort!` accepts.
- `Pkg.test()` then errored with `Package Random not found in current path` because the test-only `Random` import wasn't declared in the `[extras]` / `[targets]` table. Added it alongside `Test`.

Both failures were straightforward and resolved on the first retry. No 3-strike situations.

## Files changed

```
Project.toml                                        +4 -1   (Clustering dep, Random test-only)
Manifest.toml                                       regenerated (Clustering + transitive deps)
src/SMLMClustering.jl                               +2 -1   (export + include backends/dbscan.jl)
src/backends/dbscan.jl                              +125   (new)
test/runtests.jl                                    +3 -1   (include test_dbscan.jl; comment update)
test/test_dbscan.jl                                 +180   (new)
rounds/round_002_dbscan-backend.md                  new
STATUS.md                                           Current State + Priority 2 → DONE + Round History
KNOWLEDGE_BASE.md                                   V5 (kwdef config pattern) + V6 (eps_nm/μm conversion rule)
```

## Confidence

High. The blob-correctness test is a direct check on DBSCAN's behavior with known ground truth (σ=10 nm clusters at 1-μm separation under eps=100 nm is firmly in the high-signal regime; a spurious extra cluster or merge would be an algorithm bug). Per-dataset namespacing, `remove_unclustered`, 3D dispatch, and error paths each have dedicated tests, so regressions from subsequent backend rounds will surface immediately.

One item to watch in later rounds: `cluster_sizes` under `per_dataset=true` is a flat concatenation of per-dataset size vectors — index-by-id only holds within a dataset's range, not across the full array. Tests assert the sum invariant but don't lock in a specific per-index layout. If any consumer (e.g. SMLMAnalysis's `analyze()` wrapper) wants dataset-aware size lookup, we'll need a small structural change; noted as a future consideration but not a blocker for HDBSCAN.

## External consultations (if any)

None.

## Next steps

Priority 3: `HDBSCANConfig`. Same shared-fields shape (V3), add `min_cluster_size::Int` as algorithm-specific (HDBSCAN's `min_cluster_size` is semantically distinct from DBSCAN's `min_points`-doubled role — HDBSCAN separates "min_samples" from "min_cluster_size"). Backend library: Clustering.jl provides `hdbscan`; need to verify its public surface and how it reports noise. The test shape can mirror `test_dbscan.jl` closely: synthetic-blob correctness, label writing, per-dataset namespacing, remove_unclustered, 3D path, argument validation, empty SMLD.

## Questions posted

None.
