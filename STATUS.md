# STATUS — SMLMClustering

*Live project state. Read in Phase 1 of every round. Updated in Phase 4.*

---

## Current State

All three active backends (`DBSCANConfig`, `HierarchicalConfig`, `VoronoiConfig`) are working end-to-end with full documentation. `cluster(smld, cfg)` is now **non-mutating** — input emitters are deep-copied, labels are written onto the copy (KB V9). `HierarchicalConfig` renamed its cut-height field to `cut_threshold` (linkage-dependent unit: nm for distance linkages, μm²-ish for Ward) and added optional `n_clusters::Union{Int,Nothing}=nothing`; exactly one must be set (KB V6 caveat + V9). `ClusterInfo.cluster_sizes` stays as a flat `Vector{Int}` per Q4. `README.md` and `api_overview.md` updated to match. Full test suite is **175/175 passing in 29.1 s**. Q3/Q4/Q5 are all PROCESSED. The one remaining external blocker is Priority 7 (Keith creating the `JuliaSMLM/SMLMClustering.jl` GitHub org repo so @analysis can wire in via `[sources]`); otherwise the package is feature-complete for the 1.0 lineup. Priorities 11/12/13/14 (`cluster_statistics` + `HopkinsConfig`) are now unblocked for autonomous rounds.

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
11. [HIGH] Scaffold the `cluster_statistics` sibling interface: `abstract type AbstractStatisticsConfig <: SMLMData.AbstractSMLMConfig`, `struct ClusterStatisticsInfo <: SMLMData.AbstractSMLMInfo` with fields `n_locs_in::Int`, `statistic::Float64`, `statistic_name::Symbol`, `algorithm::Symbol`, `elapsed_s::Float64`, and `extras::Dict{Symbol,Any}` (for backends that produce multi-value stats — most fill just `statistic`), and a fallback `cluster_statistics(smld::BasicSMLD, cfg::AbstractStatisticsConfig)` that errors pointing at concrete backends. Return shape: `(smld, ClusterStatisticsInfo)` — SMLD is the SAME reference as the input (passthrough, no copy); this semantic is intentional and different from `cluster()`'s copy-semantics guarantee. Docstring must call this out explicitly so users don't assume symmetry across the two entry points. **Convention for backends that produce vector-valued outputs** (e.g. per-cluster silhouette scores, cluster-size distribution): put a summary scalar in `statistic` (overall mean silhouette, median cluster size, etc.) and the full vector in `extras` under a descriptive key. Document this convention in the `cluster_statistics` docstring so backend authors don't reinvent it. Export `AbstractStatisticsConfig`, `ClusterStatisticsInfo`, `cluster_statistics`. Tests in `test/runtests.jl` mirror the Round 001 pattern (subtyping claims, fallback error). Record the new abstract hierarchy + entry point + passthrough-same-reference + summary-scalar-plus-extras convention as KB V7. — TODO
12. [HIGH] Implement `HopkinsConfig <: AbstractStatisticsConfig` — computes the Hopkins statistic (clustering tendency) on the unlabeled coord set. Shared-ish fields: `use_3d::Bool=false`, `per_dataset::Bool=true` (per-dataset computation with a reported aggregate or a vector — design sub-question). Algorithm-specific: `n_samples::Int=20` (or `0.05` fraction), `seed::Union{Int,Nothing}=nothing` for reproducibility, `random_repeats::Int=1` (average over repeats to reduce variance). `cluster_statistics(smld, ::HopkinsConfig)` fills `ClusterStatisticsInfo(algorithm=:hopkins, statistic_name=:hopkins, statistic=H)`. Uses existing `_coords_matrix` helper; samples need uniform reference bounding box from the coord extrema. — TODO
13. [MEDIUM] Test suite for `cluster_statistics`: synthetic uniform-random SMLD → Hopkins ≈ 0.5; synthetic tight-blob SMLD → Hopkins close to 1.0; per-dataset handling; argument validation (reject labeled input? probably accept it and ignore labels since `per_dataset` uses `dataset` field); empty SMLD. Mirror the `test_dbscan.jl` structure. — TODO
14. [LOW] Update README + `api_overview.md` to cover the `cluster_statistics` entry point alongside `cluster`, explain the split (labeling vs diagnostics), and document `HopkinsConfig`. — TODO

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
| 009 | Idle close + surface @analysis Q3-Q5 input + add Stop Condition | opus | done | No unblocked priorities; @analysis consumer perspective inlined on Q3/Q4/Q5; Stop Condition halts autonomous rounds until Keith unblocks |
| 010 | Stop Condition halt | opus | stopped | Phase 2 halted on the active Stop Condition; no source changes; round-history stamp only |
| 011 | Process ANSWERED Q3/Q4/Q5 | opus | done | Non-mutating cluster (deepcopy) + HierarchicalConfig cut_threshold/n_clusters rename & split; Q4 no-op; KB V6 caveat + V9; 175/175 tests pass |

---

## Stop Conditions

<!-- Halt Phase 2 before any work if any of these are met. Add a line to halt autonomous operation without touching the lock. -->

- (none — Q3/Q4/Q5 processed in Round 011; Priorities 11–14 are now unblocked for autonomous work. Priority 7 remains blocked on Keith creating the JuliaSMLM org repo, but does not gate all rounds since other priorities can progress.)

<!-- Examples:
- `Test suite broken on main branch (block all rounds until fixed)`
- `Paused by kalidke on 2026-04-15: waiting for external review`
- `All priorities DONE — project complete`
-->
