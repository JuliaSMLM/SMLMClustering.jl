# STATUS — SMLMClustering

*Live project state. Read in Phase 1 of every round. Updated in Phase 4.*

---

## Current State

<!-- One paragraph: what works end-to-end right now, what is in flight, what recent rounds established. -->
<!-- Example: "Scalar PSF training validated at 98.4% accuracy on R044. Vector PSF model partially wired — forward pass works, training loop untested. Reactant backend default; Enzyme disabled pending upstream fix." -->

Interface + first backend are working. `Project.toml` depends on `SMLMData v0.7` and `Clustering v0.15`. `src/types.jl` defines `AbstractClusterConfig <: SMLMData.AbstractSMLMConfig` and `ClusterInfo <: SMLMData.AbstractSMLMInfo`. `src/backends/dbscan.jl` defines `DBSCANConfig(eps_nm, min_points=5, use_3d=false, per_dataset=true, remove_unclustered=false)` and a specialized `cluster(smld, ::DBSCANConfig)` method that converts nm→μm, groups by dataset when requested, runs `Clustering.dbscan`, writes labels to `emitter.id` (0=noise, 1..K=cluster with per-dataset-local namespaces), and returns `(smld_out, ClusterInfo(:dbscan, ...))`. Full suite is 56/56 passing in 2.2 s and covers three-blob correctness, per-dataset namespace locality, `remove_unclustered` filtering, 3D dispatch, argument validation, and empty SMLDs. SMLMAnalysis's `const DBSCANConfig = SMLMClustering.DBSCANConfig` re-export is now safe. Next: HDBSCAN (Priority 3).

---

## Active Threads

<!-- Work explicitly in flight across multiple rounds. Each thread is a line: name, current state, owning round. Keep short. -->

(none)

---

## Future Priorities

<!-- Ordered list. Severity tags [CRITICAL] | [HIGH] | [MEDIUM] | [LOW] optional. Status tags TODO | IN PROGRESS | BLOCKED | DONE. -->
<!-- Human seeds these on init. Rounds reorder / update status / add discovered items, never silently drop. -->

1. [HIGH] Scaffold package skeleton: add SMLMData dep to Project.toml, define `AbstractClusterConfig <: SMLMData.AbstractSMLMConfig`, define `ClusterInfo <: AbstractSMLMInfo` with fields `n_locs_in, n_clustered, n_noise, n_clusters, cluster_sizes::Vector{Int}, algorithm::Symbol, elapsed_s`, define the `cluster(smld::BasicSMLD, cfg::AbstractClusterConfig) -> (smld, ClusterInfo)` entry point with no backends wired yet — DONE (Round 001)
2. [HIGH] Implement `DBSCANConfig` backend: shared fields (`min_points`, `use_3d`, `per_dataset`, `remove_unclustered`) + algorithm field (`eps_nm`); dispatch `cluster(smld, ::DBSCANConfig)` writing labels to `emitter.id` with per-dataset handling — DONE (Round 002)
3. [HIGH] Implement `HDBSCANConfig` backend: same shared fields + `min_cluster_size` — TODO
4. [MEDIUM] Implement `VoronoiConfig` backend: density-based via Voronoi tessellation; shared fields + `density_factor` — TODO
5. [MEDIUM] Implement `HierarchicalConfig` backend: shared fields + `cut_nm` threshold — TODO
6. [MEDIUM] Test suite: one test file per backend covering clustering correctness on a known-label synthetic SMLD, per-dataset handling, and `remove_unclustered` behavior — IN PROGRESS (DBSCAN covered in Round 002; HDBSCAN/Voronoi/Hierarchical pending)
7. [LOW] API overview + README covering the four backends and the `(smld, ClusterInfo)` tuple convention — TODO

---

## Round History

| NNN | Focus | Status | Key finding |
|-----|-------|--------|-------------|
| 000 | Initial scaffold | done | Round system installed via /round-init; scope and interface agreed with @analysis |
| 001 | Scaffold package skeleton | done | SMLMData dep + AbstractClusterConfig + ClusterInfo + cluster() fallback landed; 11/11 tests pass |
| 002 | DBSCAN backend | done | DBSCANConfig + cluster dispatch via Clustering.dbscan; per-dataset namespacing verified; 56/56 tests pass |

---

## Stop Conditions

<!-- Halt Phase 2 before any work if any of these are met. Add a line to halt autonomous operation without touching the lock. -->

- (none)

<!-- Examples:
- `Test suite broken on main branch (block all rounds until fixed)`
- `Paused by kalidke on 2026-04-15: waiting for external review`
- `All priorities DONE — project complete`
-->
