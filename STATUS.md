# STATUS — SMLMClustering

*Live project state. Read in Phase 1 of every round. Updated in Phase 4.*

---

## Current State

All three active backends (`DBSCANConfig`, `HierarchicalConfig`, `VoronoiConfig`) are working end-to-end with full documentation. `README.md` covers the entry point, all three backends with examples, the `ClusterInfo` field table, and shared config fields. `api_overview.md` provides LLM-parseable API reference for all public exports. Full test suite is **159/159 passing in 28.6 s**. Three OPEN questions remain for Keith: the `cluster` vs `cluster!` naming convention (Q3), whether `ClusterInfo.cluster_sizes` should track dataset provenance (Q4), and the Ward-linkage + `cut_nm` unit mismatch on `HierarchicalConfig` (Q5). The package is feature-complete for the 1.0 lineup pending the GitHub push (Priority 7, blocked on Keith creating the org repo) and resolution of Q3/Q4/Q5.

---

## Active Threads

<!-- Work explicitly in flight across multiple rounds. Each thread is a line: name, current state, owning round. Keep short. -->

(none)

---

## Future Priorities

<!-- Ordered list. Severity tags [CRITICAL] | [HIGH] | [MEDIUM] | [LOW] optional. Status tags TODO | IN PROGRESS | BLOCKED | DONE. -->
<!-- Human seeds these on init. Rounds reorder / update status / add discovered items, never silently drop. -->

1. [HIGH] Scaffold package skeleton — DONE (Round 001)
2. [HIGH] Implement `DBSCANConfig` backend — DONE (Round 002)
3. [HIGH] Implement `HDBSCANConfig` backend: same shared fields + `min_cluster_size` — BLOCKED (no lightweight Julia HDBSCAN library; see D1 in KNOWLEDGE_BASE.md). Q1 resolved: keep the slot open and revisit when a suitable library appears.
4. [MEDIUM] Implement `VoronoiConfig` backend: density-based via Voronoi tessellation; shared fields + `density_factor` — DONE (Round 004)
5. [MEDIUM] Implement `HierarchicalConfig` backend: shared fields + `cut_nm` threshold — DONE (Round 003)
6. [MEDIUM] Test suite: one test file per backend covering clustering correctness on a known-label synthetic SMLD, per-dataset handling, and `remove_unclustered` behavior — IN PROGRESS (DBSCAN + Hierarchical + Voronoi covered; HDBSCAN pending a library)
7. [HIGH] Push SMLMClustering to the `JuliaSMLM` GitHub org so SMLMAnalysis can pull it in via `[sources]` with `rev="main"`. @analysis's Q2 answer confirms the integration mode; @analysis is waiting on the URL + branch name. Needs Keith to create the repo under the org (external action) — TODO
8. [LOW] API overview + README covering the four backends and the `(smld, ClusterInfo)` tuple convention — DONE (Round 008)
9. [MEDIUM] Replace hand-written `_pairwise_distances` in `src/utils.jl` with `Distances.pairwise(Euclidean(), X; dims=2)`. `Distances.jl` is already transitively available via Clustering.jl; adding it to `[deps]` is cheap. Gives BLAS-backed pairwise for the hierarchical backend on large groups. Source: Round 005 review I4 (AGREE-substantial). — DONE (Round 006)
10. [MEDIUM] Guard Voronoi backend against duplicate-coordinate inputs. `DelaunayTriangulation.get_area` raises `KeyError` on exact-coincident generators; the current Voronoi path will crash rather than error cleanly. Either deduplicate before `triangulate` or wrap with `ArgumentError`. Add a regression test at `test/test_voronoi.jl` for duplicate coords. Source: Round 005 review I6 (AGREE-substantial). — DONE (Round 007)

---

## Round History

| NNN | Focus | Model | Status | Key finding |
|-----|-------|-------|--------|-------------|
| 000 | Initial scaffold | opus | done | Round system installed via /round-init; scope and interface agreed with @analysis |
| 001 | Scaffold package skeleton | opus | done | SMLMData dep + AbstractClusterConfig + ClusterInfo + cluster() fallback landed; 11/11 tests pass |
| 002 | DBSCAN backend | opus | done | DBSCANConfig + cluster dispatch via Clustering.dbscan; per-dataset namespacing verified; 56/56 tests pass |
| 003 | Hierarchical backend | sonnet | done | HierarchicalConfig via hclust+cutree+min_points filter; HDBSCAN blocked (no library); 108/108 tests pass |
| 004 | Voronoi backend | opus | done | VoronoiConfig via DelaunayTriangulation.jl (SR-Tesseler density); 2D only; 157/157 tests pass; Q1+Q2 processed |
| 005 | Review at Round 005 | opus | done | /review-code: trivial fixes applied (helpers, time_ns, show); I4+I6 as P9+P10; Q3/Q4/Q5 posted |
| 006 | Distances.pairwise replacement | sonnet | done | _pairwise_distances → one-liner via Distances.pairwise; Distances dep explicit; 157/157 pass |
| 007 | Voronoi duplicate-coordinate guard | sonnet | done | ArgumentError before triangulate on exact-coincident (x,y) pairs; 159/159 pass |
| 008 | API overview + README | sonnet | done | README.md rewritten; api_overview.md created; 159/159 pass; package feature-complete pending GitHub push |

---

## Stop Conditions

<!-- Halt Phase 2 before any work if any of these are met. Add a line to halt autonomous operation without touching the lock. -->

- (none)

<!-- Examples:
- `Test suite broken on main branch (block all rounds until fixed)`
- `Paused by kalidke on 2026-04-15: waiting for external review`
- `All priorities DONE — project complete`
-->
