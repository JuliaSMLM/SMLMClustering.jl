# STATUS — SMLMClustering

*Live project state. Read in Phase 1 of every round. Updated in Phase 4.*

---

## Current State

All five labeling backends (`DBSCANConfig`, `HDBSCANConfig`, `HierarchicalConfig`, `VoronoiConfig`, `MRFDensityClusterConfig`) and both spatial-statistic backends (`HopkinsConfig`, `VoronoiDensityConfig`) are live end-to-end with full documentation. `cluster(smld, cfg)` is non-mutating (deep-copies emitters per KB V9); `cluster_statistics(smld, cfg)` is pass-through (returns input SMLD reference per KB V10). Task A (P11+P12+P13+P14) and Task B (P15) landed 2026-04-27 via direct execution; `MRFDensityClusterConfig` (V11 four-step pipeline: Voronoi density → multi-component GMM → multi-class Potts MRF via ICM → CC on foreground) landed 2026-04-29 to address the "missing middles + spurious smalls" pathology from genmab dSTORM data without per-dataset density tuning. **Test suite split into a fast tier (161/161, default `Pkg.test()`) and a thorough tier gated by `SMLM_TEST_FULL=true|1|yes` (981/981); shared cross-package convention with SMLMAnalysis / SMLMBaGoL / SMLMDriftCorrection.** README and api_overview cover both entry points, the labeling-vs-diagnostics split, and the summary-scalar+extras convention. Only remaining open work is Priority 7 (Keith creating the `JuliaSMLM/SMLMClustering.jl` GitHub org repo so @analysis can wire in via `[sources]`).

---

## Active Threads

<!-- Work explicitly in flight across multiple rounds. Each thread is a line: name, current state, owning round. Keep short. -->

(none — Task A + Task B both landed 2026-04-27)

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
11. [HIGH] Scaffold the `cluster_statistics` sibling interface: `abstract type AbstractStatisticsConfig <: SMLMData.AbstractSMLMConfig`, `struct ClusterStatisticsInfo <: SMLMData.AbstractSMLMInfo` with fields `n_locs_in::Int`, `statistic::Float64`, `statistic_name::Symbol`, `algorithm::Symbol`, `elapsed_s::Float64`, and `extras::Dict{Symbol,Any}` (for backends that produce multi-value stats — most fill just `statistic`), and a fallback `cluster_statistics(smld::BasicSMLD, cfg::AbstractStatisticsConfig)` that errors pointing at concrete backends. Return shape: `(smld, ClusterStatisticsInfo)` — SMLD is the SAME reference as the input (passthrough, no copy); this semantic is intentional and different from `cluster()`'s copy-semantics guarantee. Docstring must call this out explicitly so users don't assume symmetry across the two entry points. **Convention for backends that produce vector-valued outputs** (e.g. per-cluster silhouette scores, cluster-size distribution): put a summary scalar in `statistic` (overall mean silhouette, median cluster size, etc.) and the full vector in `extras` under a descriptive key. Document this convention in the `cluster_statistics` docstring so backend authors don't reinvent it. Export `AbstractStatisticsConfig`, `ClusterStatisticsInfo`, `cluster_statistics`. Tests in `test/runtests.jl` mirror the Round 001 pattern (subtyping claims, fallback error). Record the new abstract hierarchy + entry point + passthrough-same-reference + summary-scalar-plus-extras convention as KB V7. — DONE (Task A direct execution 2026-04-27; KB entry landed as V10 not V7)
12. [HIGH] Implement `HopkinsConfig <: AbstractStatisticsConfig` — computes the Hopkins statistic (clustering tendency) on the unlabeled coord set. Shared-ish fields: `use_3d::Bool=false`, `per_dataset::Bool=true` (per-dataset computation with a reported aggregate or a vector — design sub-question). Algorithm-specific: `n_samples::Int=20` (or `0.05` fraction), `seed::Union{Int,Nothing}=nothing` for reproducibility, `random_repeats::Int=1` (average over repeats to reduce variance). `cluster_statistics(smld, ::HopkinsConfig)` fills `ClusterStatisticsInfo(algorithm=:hopkins, statistic_name=:hopkins, statistic=H)`. Uses existing `_coords_matrix` helper; samples need uniform reference bounding box from the coord extrema. — DONE (Task A direct execution 2026-04-27)
13. [MEDIUM] Test suite for `cluster_statistics`: synthetic uniform-random SMLD → Hopkins ≈ 0.5; synthetic tight-blob SMLD → Hopkins close to 1.0; per-dataset handling; argument validation (reject labeled input? probably accept it and ignore labels since `per_dataset` uses `dataset` field); empty SMLD. Mirror the `test_dbscan.jl` structure. — DONE (Task A direct execution 2026-04-27)
14. [LOW] Update README + `api_overview.md` to cover the `cluster_statistics` entry point alongside `cluster`, explain the split (labeling vs diagnostics), and document `HopkinsConfig`. — DONE (Task A direct execution 2026-04-27)
15. [MEDIUM] Per-emitter Voronoi density utility (e.g. `voronoi_density(smld; use_3d=false, per_dataset=true) -> Vector{Float64}` returning per-emitter ρ_i = 1/A_i, μm⁻²). Currently the area computation already lives inside `VoronoiConfig`'s backend (src/backends/voronoi.jl ~line 121) but is consumed for the density threshold and discarded — exposing it as a sibling utility lets downstream callers run their own thresholding (Otsu, GMM on log ρ_i, fixed cutoff, etc.) for cell-structure masking. Likely slots into the `cluster_statistics` interface (P11) as `VoronoiDensityConfig` with the per-emitter vector in `extras[:density_per_emitter]` and a summary scalar (median density?) in `statistic`. Mirror voronoi.jl's degeneracy guards (groups <3, duplicate coords). Motivation: paper-genmab-hexabody A431 cell-structure masking pre-step (consumer is project-local against DelaunayTriangulation.jl in the interim). — DONE (Task B direct execution 2026-04-27; landed as `VoronoiDensityConfig <: AbstractStatisticsConfig`, `extras[:density_per_emitter]` + `extras[:area_per_emitter]` flat in original emitter order, summary `statistic = median(non-NaN densities)`)
16. [LOW] Cross-package adoption of the test-tier split convention. SMLMClustering shipped the convention 2026-04-27 (env-var `SMLM_TEST_FULL`, permissive truthy-check accepting `true|1|yes` case-insensitively; default-off thorough tier; `@info` skip-message when off; per-testset `if SMLM_TEST_FULL ... end` gating in each test_*.jl file). Pattern is research-validated against Flux.jl precedent (`<PKG>_TEST_<DOMAIN>` shape) and ratified by @analysis. Shared name `SMLM_TEST_FULL` (no per-package namespacing) is intentional — applies across the SMLM family. **Pending propagation to: SMLMAnalysis (@analysis primed and willing — currently 98 tests, ~4s, no thorough tier yet so the split mostly carves out a place to put thorough tests), SMLMBaGoL (agent dead, needs revival), SMLMDriftCorrection (agent dead, needs revival).** — TODO
17. [LOW] v2 enhancements for `MRFDensityClusterConfig`: (a) graph-cuts MAP inference (binary or α-expansion) replacing ICM as default — `inference=:graph_cuts` symbol slot is reserved in v1 but raises `ArgumentError`; would require either a new dep (BoykovKolmogorov.jl, GraphsFlows.jl) or an in-house max-flow implementation. (b) Soft-posterior output mode — store `extras[:posterior_per_emitter]::Matrix{Float64}` (n × n_regimes) when a flag is set; useful for downstream uncertainty-aware aggregations. (c) k-NN density estimator alternative to Voronoi — `density_estimator=:knn` with a `density_k` parameter; useful when Voronoi tessellation is undesirable (e.g. very large groups, ill-conditioned hull). All three would be additive — current pipeline stays the default. — TODO

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
| Task A — direct | cluster_statistics + HopkinsConfig (P11+P12+P13+P14 bundled) | opus | done | Sibling pass-through interface + Hopkins backend (KDTree NN via NearestNeighbors.jl) + tests + docs + KB V10; bypassed round dispatcher per Keith's hour-scale push for paper-genmab-hexabody v1; 221/221 tests pass |
| Task B — direct | VoronoiDensityConfig (P15) | opus | done | Per-emitter Voronoi density via DelaunayTriangulation; flat per-emitter density+area vectors in extras, median density as summary; mirrors voronoi.jl degeneracy guards; bypassed round dispatcher (same paper push as Task A); follows V10 pass-through + summary-scalar+extras conventions; 861/861 tests pass |
| Test-tier split — direct | SMLM_TEST_FULL convention (P16) | opus | done | Cross-package convention adopted: `SMLM_TEST_FULL` env var (permissive `true\|1\|yes` truthy-check) gates thorough tier; default-off; `@info` skip-message; per-testset `if SMLM_TEST_FULL ... end` gating in every test_*.jl. Fast tier 101/101 in 30s, thorough tier 861/861 in 31s; both `SMLM_TEST_FULL=true` and `SMLM_TEST_FULL=1` verified equivalent. Research-validated against Flux.jl precedent; shape ratified by @analysis. Pending propagation to SMLMAnalysis / SMLMBaGoL / SMLMDriftCorrection. |
| MRF density backend — direct | MRFDensityClusterConfig (V11) | opus | done | Adaptive-density clustering pipeline: Voronoi density → n-component 1D GMM (sorted ascending by mean, regime 1 = lowest = noise) → multi-class Potts MRF via ICM with auto-λ (MAD of unary range) → BFS connected-components on foreground (regime ≥ 2) over Delaunay (default) or kNN graph. `regime_thresholds` override bypasses GMM. Per-cluster outputs in `emitter.id`; metadata stamps `mrf_regime_per_emitter` / `mrf_lambda_used` / `mrf_regime_means`. Refactored `_voronoi_areas` into utils.jl (shared with VoronoiDensityConfig). Fast tier 161/161 in 33s; thorough tier 981/981 in 34.7s. KB V11. P17 logged for v2 enhancements (graph-cuts, soft posteriors, kNN density estimator). |

---

## Stop Conditions

<!-- Halt Phase 2 before any work if any of these are met. Add a line to halt autonomous operation without touching the lock. -->

- (none — Q3/Q4/Q5 processed in Round 011; Priorities 11–14 are now unblocked for autonomous work. Priority 7 remains blocked on Keith creating the JuliaSMLM org repo, but does not gate all rounds since other priorities can progress.)

<!-- Examples:
- `Test suite broken on main branch (block all rounds until fixed)`
- `Paused by kalidke on 2026-04-15: waiting for external review`
- `All priorities DONE — project complete`
-->
