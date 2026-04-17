# Round 003 — Hierarchical backend

**Date:** 2026-04-17
**Status:** done
**Priority worked on:** Priority 5 — Implement `HierarchicalConfig` backend (worked on instead of Priority 3 due to HDBSCAN library gap; see Dead End D1)

## Hypothesis

With the DBSCAN backend pattern established, adding a hierarchical backend using `Clustering.hclust` + `Clustering.cutree` should be straightforward: build an O(n²) distance matrix, cut the dendrogram at `cut_nm`, filter small clusters as noise via `min_points`, and write labels to `emitter.id`. The existing `_coords_matrix` helper can be factored into a shared utils.jl to avoid duplication.

## What was attempted

- **HDBSCAN library research**: checked `Clustering.jl` 0.15.8 — no `hdbscan` symbol. Searched Julia General registry — no registered HDBSCAN package. Found `baggepinnen/HDBSCAN.jl` (GitHub-only, PyCall wrapper requiring Python hdbscan) and `HorseML.jl` (registered, pure Julia HDBSCAN but pulls in CUDA + NNlib + Zygote — too heavy for a pure algorithm package per V4). Priority 3 declared BLOCKED; D1 Dead End documented.

- **Pivoted to Priority 5 (HierarchicalConfig)**: `Clustering.hclust` + `Clustering.cutree` already in our dependency. Tested edge cases (1-point, 2-point, empty). Confirmed that cutting at 100 nm on σ=10 nm blobs + min_points=5 filtering correctly recovers 3 clusters from 33 raw hierarchical clusters.

- **`src/utils.jl`** (new): extracted `_coords_matrix` and added `_pairwise_distances`. Both shared by DBSCAN and Hierarchical backends. Loaded before backend files via `SMLMClustering.jl`.

- **`src/backends/dbscan.jl`**: removed the private `_coords_matrix` definition (now in utils.jl); added a comment referencing utils.jl. All 45 DBSCAN tests continue to pass.

- **`src/backends/hierarchical.jl`** (new): `HierarchicalConfig` with fields `cut_nm::Float64`, `linkage::Symbol=:ward`, plus shared fields (V5/V3). `cluster` method: groups by dataset, builds pairwise Euclidean distance matrix (μm), calls `hclust` + `cutree(h=cut_nm/1000)`, applies min_points filter to relabel small clusters as noise, writes local label namespace per V3.

- **`src/SMLMClustering.jl`**: added `include("utils.jl")` and `include("backends/hierarchical.jl")`; exported `HierarchicalConfig`.

- **`test/test_hierarchical.jl`** (new, 52 tests): mirrors `test_dbscan.jl` structure. Testsets: config construction, three-blob correctness (single linkage), labels + remove_unclustered, per_dataset namespace, argument validation (cut_nm=0, min_points=0, invalid linkage), use_3d on 2D error, 3D path, empty SMLD, min_points filtering as noise.

- **`test/runtests.jl`**: added `include("test_hierarchical.jl")`.

- **First test run**: 107/108 passed. One failure: per_dataset test showed dataset 2 got ids [3,4] instead of [1,2]. Root cause: `label_map[orig] = length(cluster_sizes)` accumulated globally across groups instead of resetting per group. Fixed by adding a `k_local` counter that resets each iteration.

## What worked

- Full suite green: **108/108 passing in 4.8 s**.
- DBSCAN's 45 tests continue to pass after the `_coords_matrix` refactor to utils.jl.
- Three-blob test with single linkage at 100 nm correctly recovers `n_clusters=3` from 33 raw clusters; 30 scattered noise points become singletons filtered to noise by `min_points=5`.
- Per-dataset label locality: two datasets each with two blobs produce ids `[1, 2]` within each dataset (4 total), not `[1, 2, 3, 4]`.
- The `min_points`-as-noise mechanism works: blobs of 10 pts survive at min_points=5 but become noise at min_points=15.
- 3D path tested (two z-separated blobs cluster correctly under use_3d=true).

## What failed

- **First label-namespace bug**: global `length(cluster_sizes)` used as label counter across per-dataset groups, producing non-local ids. Fixed in one retry by adding a per-group `k_local` counter.

## Files changed

```
src/utils.jl                                +47    (new — _coords_matrix + _pairwise_distances)
src/backends/dbscan.jl                      -26 +5  (remove _coords_matrix; add utils.jl comment)
src/backends/hierarchical.jl               +117   (new)
src/SMLMClustering.jl                       +3 -2  (include utils + hierarchical; export HierarchicalConfig)
test/test_hierarchical.jl                  +200   (new — 52 tests)
test/runtests.jl                            +1     (include test_hierarchical.jl)
rounds/round_003_hierarchical-backend.md   new
STATUS.md                                  updated
KNOWLEDGE_BASE.md                          +D1
```

## Confidence

High. 108/108 tests pass including 3D path, per-dataset namespace, remove_unclustered, min_points filter, and argument validation. The per-dataset namespace bug was caught immediately by the test and fixed in one retry.

One forward note: the `cluster_sizes` vector for hierarchical clustering with `per_dataset=true` is a flat concatenation of per-dataset size vectors (same caveat as DBSCAN Round 002 — index-by-id is only valid within a dataset's range). If SMLMAnalysis needs dataset-aware size lookup, a structural change to `ClusterInfo` would be needed.

## External consultations (if any)

- WebSearch: "HDBSCAN.jl Julia package GitHub 2024 2025" — found `baggepinnen/HDBSCAN.jl` (PyCall wrapper, needs Python), `HorseML.jl` (pure Julia, too heavy). Used to rule out Priority 3 and justify pivot to Priority 5.

## Next steps

Priority 3 (HDBSCAN) is BLOCKED pending a lightweight library. Priority 5 is now DONE. Next unblocked priority is Priority 4 (VoronoiConfig, MEDIUM) — Voronoi tessellation backend. This requires a Voronoi library; needs investigation in Phase 1 of the next round.

## Questions posted

None.
