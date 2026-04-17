# STATUS — SMLMClustering

*Live project state. Read in Phase 1 of every round. Updated in Phase 4.*

---

## Current State

Two backends are working end-to-end. `src/utils.jl` provides shared helpers (`_coords_matrix`, `_pairwise_distances`). `src/backends/dbscan.jl` implements `DBSCANConfig` (eps_nm, Clustering.dbscan). `src/backends/hierarchical.jl` implements `HierarchicalConfig` (cut_nm, linkage, agglomerative via `Clustering.hclust` + `cutree`, with `min_points` noise filtering). Full suite is 108/108 passing in 4.8 s. HDBSCAN (Priority 3) is BLOCKED: `Clustering.jl` 0.15.8 has no `hdbscan`, the only registered Julia HDBSCAN is in `HorseML.jl` (CUDA + NNlib + Zygote — too heavy), and the GitHub-only `baggepinnen/HDBSCAN.jl` requires Python. Next viable priority is Priority 4 (VoronoiConfig, MEDIUM).

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
3. [HIGH] Implement `HDBSCANConfig` backend: same shared fields + `min_cluster_size` — BLOCKED (no lightweight Julia HDBSCAN library; see D1 in KNOWLEDGE_BASE.md)
4. [MEDIUM] Implement `VoronoiConfig` backend: density-based via Voronoi tessellation; shared fields + `density_factor` — TODO
5. [MEDIUM] Implement `HierarchicalConfig` backend: shared fields + `cut_nm` threshold — DONE (Round 003)
6. [MEDIUM] Test suite: one test file per backend covering clustering correctness on a known-label synthetic SMLD, per-dataset handling, and `remove_unclustered` behavior — IN PROGRESS (DBSCAN + Hierarchical covered; Voronoi/HDBSCAN pending)
7. [LOW] API overview + README covering the four backends and the `(smld, ClusterInfo)` tuple convention — TODO

---

## Round History

| NNN | Focus | Model | Status | Key finding |
|-----|-------|-------|--------|-------------|
| 000 | Initial scaffold | opus | done | Round system installed via /round-init; scope and interface agreed with @analysis |
| 001 | Scaffold package skeleton | opus | done | SMLMData dep + AbstractClusterConfig + ClusterInfo + cluster() fallback landed; 11/11 tests pass |
| 002 | DBSCAN backend | opus | done | DBSCANConfig + cluster dispatch via Clustering.dbscan; per-dataset namespacing verified; 56/56 tests pass |
| 003 | Hierarchical backend | sonnet | done | HierarchicalConfig via hclust+cutree+min_points filter; HDBSCAN blocked (no library); 108/108 tests pass |

---

## Stop Conditions

<!-- Halt Phase 2 before any work if any of these are met. Add a line to halt autonomous operation without touching the lock. -->

- (none)

<!-- Examples:
- `Test suite broken on main branch (block all rounds until fixed)`
- `Paused by kalidke on 2026-04-15: waiting for external review`
- `All priorities DONE — project complete`
-->
