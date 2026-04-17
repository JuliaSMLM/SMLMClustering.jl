# Round 005 — Review at Round 005

**Date:** 2026-04-17
**Status:** done
**Priority worked on:** Forced review cycle (review_cycle_rounds = 5). Normal priority work replaced by `/review-code` + classification per `.claude/round/review-protocol.md`.

## Hypothesis

With three backends shipped, the first 5-round review should surface (a) duplication across the backends that accumulated as each one was copy-pasted from the previous, (b) any brittle corners from the rapid implementation (duplicate-coord Voronoi, `time()` monotonicity, Ward+`cut_nm` unit mismatch), and (c) interface-level questions now that the four-backend surface is near-final. The review should find mostly trivial/substantial findings; the deeper design questions get deferred to QUESTIONS.md for Keith.

## What was attempted

- **Invoked `/review-code`** as a comprehensive in-tree review covering file hierarchy, type hierarchy, workflow, code quality, dependencies, idiomatic Julia improvements, and strengths. Output enumerated nine findings (I1–I6, R1–R2, M1).

- **Classified each finding** per `review-protocol.md` into AGREE-trivial / AGREE-substantial / DISAGREE / DEFER:

  | Finding | Bucket | Action taken |
  |---|---|---|
  | I1. Dead `using Clustering` in `src/utils.jl` | AGREE-trivial | Fixed |
  | I2. Duplicated group-by-dataset / compact-relabel / output-SMLD blocks | AGREE-trivial | Fixed |
  | I3. `time()` → `time_ns()` for monotonic `elapsed_s` | AGREE-trivial | Fixed |
  | I4. `Distances.pairwise` replaces hand loop in `_pairwise_distances` | AGREE-substantial | New Priority 9 |
  | I5. `Base.show(::IO, ::ClusterInfo)` | AGREE-trivial | Fixed |
  | I6. Voronoi duplicate-coordinate guard (`KeyError` latent crash) | AGREE-substantial | New Priority 10 |
  | R1. `cluster` → `cluster!` rename | DEFER | Q3 OPEN |
  | R2. `ClusterInfo.cluster_sizes` dataset provenance | DEFER | Q4 OPEN |
  | M1. Ward linkage + `cut_nm` foot-gun | DEFER | Q5 OPEN |

- **Applied the trivial fixes in order:**
  1. Added three helpers to `src/utils.jl`: `_group_by_dataset(smld, per_dataset)`, `_compact_relabel!(cluster_sizes, raw_counts, min_points)` (returns `(label_map, n_added)`), and `_build_output(smld, remove_unclustered)`.
  2. Refactored all three backends to call the helpers. Removed ~40 lines of duplication.
  3. Switched all three backends' timing to `time_ns()` with a `/1e9` to yield `elapsed_s::Float64`.
  4. Added a custom `Base.show(::IO, ::ClusterInfo)` to `src/types.jl` producing `ClusterInfo(85/100 clustered, 4 clusters, algorithm=:dbscan, 123.0 ms)`.

- **Caught a test-suite regression mid-round and fixed it:** the `using Clustering` I removed from `utils.jl` had been the only module-level `using` for Clustering — the backend files reference `Clustering.dbscan` / `Clustering.hclust` / `Clustering.cutree` unqualified. Moved `using Clustering` to the module file `src/SMLMClustering.jl` where it belongs. Single strike, not a Dead End.

- **Posted Q3 (cluster vs cluster!), Q4 (cluster_sizes dataset-awareness), Q5 (Ward + cut_nm)** as OPEN items in `QUESTIONS.md`. Each one follows the four-block shape with a Short Question answerable by Keith without reading this round file, and a Technical detail block for the round that processes the answer.

- **Updated STATUS.md:** replaced `## Current State` with a fresh summary of the new reality; appended Round 005 row to Round History with Model=opus; added Priorities 9 (Distances.pairwise) and 10 (Voronoi duplicate-coord guard).

## What worked

- **Test suite fully green after refactor:** 157/157 passing in 28.5 s, identical count to Round 004. No semantic regression.
- **Code footprint shrank:** three backend files are each ~10 lines shorter after the helper extraction; utils.jl gained ~45 lines net but all three call sites became mechanical.
- **The three helpers caught all three duplicated patterns cleanly** — no awkward parameter threading required. `_compact_relabel!` returning `(label_map, added)` lets the caller accumulate `n_clustered` without an extra pass.
- **`time_ns()` migration was mechanical** — three backend call sites, one pattern change. No observable behavior difference in tests.

## What failed

- **Test-suite regression after removing `using Clustering` from utils.jl** (one strike): the backends relied on the transitive import via `utils.jl` being included first. Fixed by moving the import to the module file `src/SMLMClustering.jl`. Correct resolution — the module file is the right place for a package-wide `using` statement, not a utility file that doesn't use the symbol itself. Not a Dead End because the fix was immediate and consistent with Julia convention.

## Files changed

```
src/SMLMClustering.jl                     +1       (using Clustering at module level)
src/types.jl                              +7       (Base.show for ClusterInfo)
src/utils.jl                              +47 -2   (3 new helpers; removed dead using Clustering)
src/backends/dbscan.jl                    +4 -25   (use helpers; time_ns; _build_output)
src/backends/hierarchical.jl              +4 -33   (use helpers; time_ns; _compact_relabel!; _build_output)
src/backends/voronoi.jl                   +3 -28   (use helpers; time_ns; _compact_relabel!; _build_output)
STATUS.md                                 Current State + P9 + P10 + Round History row
QUESTIONS.md                              +Q3 (cluster vs cluster!) +Q4 (cluster_sizes) +Q5 (Ward+cut_nm)
rounds/round_005_review-at-round-005.md   new
```

## Confidence

High. The refactor preserved all 157 tests; no edge case changed. `time_ns()` is strictly better than `time()` for interval measurement. The three helpers are mechanical extractions of code that was already working — no new logic introduced. The regression during the round (removing `using Clustering` broke the backends) was detected by the test suite and fixed cleanly in one iteration.

Caveat: the `Base.show` output for `ClusterInfo` is not covered by a test. Since the struct itself is unchanged (show is pure formatting), this is acceptable; a test would be stylistic.

## External consultations (if any)

None. Review findings were all in-tree judgment calls informed by reading the backend code, tests, and KB entries.

## Next steps

- **Priority 9 (MEDIUM): Distances.pairwise replacement.** One-file edit in `utils.jl` plus `Distances` dep add. Mechanical enough to be a good Sonnet round.
- **Priority 10 (MEDIUM): Voronoi duplicate-coord guard.** One-file edit in `voronoi.jl` + one test case. Also mechanical.
- **Priority 7 (HIGH, still blocked on Keith):** GitHub push. No change this round.
- **Three OPEN questions** (Q3/Q4/Q5) await Keith's judgment before the corresponding design changes can land.

Next round should pick Priority 9 (Distances.pairwise) — top non-BLOCKED, non-Keith-gated priority. Both P9 and P10 are clean Sonnet candidates: pattern-replication against existing shape (P9 is a line swap using a library already in the transitive closure; P10 is a test + try/catch after an established pattern). Opting the next round into Sonnet via `/tmp/smlmclustering-sonnet-ok.md` below.

## Questions posted

- **Q3 — cluster vs cluster!** Should the clustering entry point be renamed to signal that it modifies its input?
- **Q4 — cluster_sizes dataset provenance** Should `ClusterInfo.cluster_sizes` carry the dataset index for each cluster?
- **Q5 — Ward + cut_nm foot-gun** Is the Ward-default + nm-cut-field combination on `HierarchicalConfig` misleading?
