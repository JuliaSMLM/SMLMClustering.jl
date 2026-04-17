# Round 006 — Distances.pairwise replacement

**Date:** 2026-04-17
**Status:** done
**Priority worked on:** Priority 9: Replace hand-written `_pairwise_distances` in `src/utils.jl` with `Distances.pairwise(Euclidean(), X; dims=2)`

## Hypothesis

`_pairwise_distances` is a 15-line hand-written nested loop computing Euclidean pairwise distances; replacing it with `Distances.pairwise` (BLAS-backed, vectorized) should be a mechanical one-liner with no semantic change and identical test outcomes. `Distances.jl` was already transitively available via `Clustering.jl` so the dep addition would be cheap.

## What was attempted

- **Added `Distances` to `Project.toml` deps** via `Pkg.add("Distances")`. Resolved immediately (already transitively present); only the explicit `[deps]` entry and a `[compat]` entry for `Distances = "0.10"` were added to Project.toml. No new packages appeared in the Manifest.
- **Added `using Distances` to `src/SMLMClustering.jl`** alongside the existing `using Clustering`.
- **Replaced `_pairwise_distances` body** with a one-liner: `Distances.pairwise(Euclidean(), X; dims=2)`. `Euclidean()` matches the hand loop's L2 norm; `dims=2` treats columns as observations, matching the `d×n` layout of `_coords_matrix`.
- **Updated stale comment in `src/backends/dbscan.jl`** that listed `_pairwise_distances` as a shared helper DBSCAN uses — DBSCAN never called it; the comment was misleading.
- **Ran full test suite:** 157/157 passing in 27.8 s, identical count to Round 005.

## What worked

- One-liner replacement is a strict drop-in: same function signature, same return type (`Matrix{Float64}` n×n symmetric), same caller interface in `hierarchical.jl:71`. No changes to hierarchical backend logic.
- Tests pass end-to-end with the new implementation — behavioral equivalence confirmed on synthetic SMLD inputs spanning 2D and 3D layouts, small and medium groups.

## What failed

Nothing. The change was as mechanical as expected.

## Files changed

```
Project.toml                              +2       (Distances dep + compat entry)
src/SMLMClustering.jl                     +1       (using Distances)
src/utils.jl                              +1 -15   (one-liner replaces 15-line loop)
src/backends/dbscan.jl                    +1 -1    (stale comment update)
STATUS.md                                 updated
rounds/round_006_distances-pairwise.md    new
```

## Confidence

High. 157/157 tests pass. The hand loop and `Distances.pairwise` are algebraically identical (L2 norm, symmetric, zero diagonal). The only caller is `hierarchical.jl`; DBSCAN was never a caller.

## External consultations (if any)

None.

## Next steps

Priority 9 DONE. Next unblocked priority is Priority 10 (MEDIUM): guard Voronoi backend against duplicate-coordinate inputs — `DelaunayTriangulation.get_area` raises `KeyError` on exact-coincident generators; need to either deduplicate before `triangulate` or wrap with `ArgumentError`, plus a regression test in `test/test_voronoi.jl`. Also a mechanical round; Sonnet-appropriate.

## Questions posted

None.
