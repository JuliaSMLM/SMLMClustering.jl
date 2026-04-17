# Round 007 — Voronoi duplicate-coordinate guard

**Date:** 2026-04-17
**Status:** done
**Priority worked on:** Priority 10: Guard Voronoi backend against duplicate-coordinate inputs

## Hypothesis

`DelaunayTriangulation.get_area` raises `KeyError` when two generators share the
same (x, y) coordinate — the triangulation silently discards the duplicate but
the per-generator area lookup by index still fires. Adding a pre-triangulation
check that throws `ArgumentError` on exact duplicates converts an opaque crash
into a descriptive error, consistent with how other degenerate inputs are handled
(`use_3d=true`, `density_factor<=0`).

## What was attempted

- **Read the Voronoi backend** (`src/backends/voronoi.jl`) to locate where `pts`
  is built and where `get_area` is called (lines 104–111). Confirmed the crash
  path: `pts` built from raw emitter coords, no deduplication, `triangulate` may
  accept duplicates silently but `get_area(vor, j)` then fails.
- **Added duplicate check** after `pts` is built (line 109 in updated file):
  `length(unique(pts)) == n || throw(ArgumentError(...))`. Single expression, no
  helper, consistent with the existing `||`-throw pattern used for `density_factor`
  and `min_points` validation.
- **Updated docstring** (`VoronoiConfig` "Degenerate input" note) to mention that
  groups with duplicate coordinates raise `ArgumentError`.
- **Added regression test** (`test/test_voronoi.jl`): two sub-cases — a single
  dataset with a duplicate point (`per_dataset=false`) and the same scenario
  with `per_dataset=true` to confirm per-group detection fires correctly.
- **Ran full test suite:** 159/159 pass in 27.7 s.

## What worked

- Check is a single line; no helper needed.
- Both sub-cases of the new testset pass: `ArgumentError` raised in the
  `per_dataset=false` group path and in the `per_dataset=true` per-group path.
- Existing 157 tests unchanged (the new cases are strictly additive).

## What failed

Nothing. The fix was as mechanical as expected.

## Files changed

```
src/backends/voronoi.jl      +8 -1   (duplicate check + updated docstring note)
test/test_voronoi.jl         +16     (new "duplicate coordinates raise ArgumentError" testset)
STATUS.md                    updated
rounds/round_007_voronoi-duplicate-guard.md  new
```

## Confidence

High. 159/159 tests pass. The new tests directly exercise the error path. The
fix fires before `triangulate` so the KeyError from `get_area` is never reached.

## External consultations (if any)

None.

## Next steps

Priority 10 DONE. Next unblocked priority is Priority 8 (LOW): API overview +
README covering the four backends and the `(smld, ClusterInfo)` tuple convention.
Sonnet-appropriate (documentation pass, no design content).

## Questions posted

None.
