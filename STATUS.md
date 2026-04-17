# STATUS — SMLMClustering

*Live project state. Read in Phase 1 of every round. Updated in Phase 4.*

---

## Current State

Three of the four planned backends are working end-to-end: `DBSCANConfig` (density-based, Clustering.jl), `HierarchicalConfig` (agglomerative, Clustering.hclust+cutree), and now `VoronoiConfig` (SR-Tesseler-style density clustering via `DelaunayTriangulation.jl`, 2D only). Shared helpers in `src/utils.jl` (`_coords_matrix`, `_pairwise_distances`) are reused across backends. Full test suite is **157/157 passing in 27.9 s** (slower than Round 003 because DelaunayTriangulation precompile + triangulate on the 180-point blob test dominates). `Project.toml` now pulls in `DelaunayTriangulation` (compat `1`) alongside `Clustering` and `SMLMData`; its transitive footprint is small (Random, ExactPredicates, AdaptivePredicates, EnumX — no plotting/CUDA/autodiff). HDBSCAN (Priority 3) remains BLOCKED per D1 — Keith's answer to Q1 was "keep the slot open, check periodically." @analysis's answer to Q2 landed late in the round: integration mode is `[sources]` with a GitHub URL on the JuliaSMLM org; this is recorded as a new priority (pushing the repo to GitHub) and is not yet done.

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
8. [LOW] API overview + README covering the four backends and the `(smld, ClusterInfo)` tuple convention — TODO

---

## Round History

| NNN | Focus | Model | Status | Key finding |
|-----|-------|-------|--------|-------------|
| 000 | Initial scaffold | opus | done | Round system installed via /round-init; scope and interface agreed with @analysis |
| 001 | Scaffold package skeleton | opus | done | SMLMData dep + AbstractClusterConfig + ClusterInfo + cluster() fallback landed; 11/11 tests pass |
| 002 | DBSCAN backend | opus | done | DBSCANConfig + cluster dispatch via Clustering.dbscan; per-dataset namespacing verified; 56/56 tests pass |
| 003 | Hierarchical backend | sonnet | done | HierarchicalConfig via hclust+cutree+min_points filter; HDBSCAN blocked (no library); 108/108 tests pass |
| 004 | Voronoi backend | opus | done | VoronoiConfig via DelaunayTriangulation.jl (SR-Tesseler density); 2D only; 157/157 tests pass; Q1+Q2 processed |

---

## Stop Conditions

<!-- Halt Phase 2 before any work if any of these are met. Add a line to halt autonomous operation without touching the lock. -->

- (none)

<!-- Examples:
- `Test suite broken on main branch (block all rounds until fixed)`
- `Paused by kalidke on 2026-04-15: waiting for external review`
- `All priorities DONE — project complete`
-->
