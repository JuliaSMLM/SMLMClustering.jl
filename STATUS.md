# STATUS — SMLMClustering

*Live project state. Read in Phase 1 of every round. Updated in Phase 4.*

---

## Current State

All five labeling backends (`DBSCANConfig`, `HDBSCANConfig`, `HierarchicalConfig`, `VoronoiConfig`, `MRFDensityClusterConfig`) and both spatial-statistic backends (`HopkinsConfig`, `VoronoiDensityConfig`) are live end-to-end with full documentation, README, and api_overview. `cluster(smld, cfg)` is non-mutating (KB V9 deepcopy); `cluster_statistics(smld, cfg)` is pass-through (KB V10 same-reference). **Round 014 (regime_thresholds calibration) closed**: new `calibrate_regime_thresholds(smld; n_regimes=2, density_estimator=:knn, density_k=20)` helper computes per-emitter log densities, fits 1D GMM, and returns analytic Bayes decision boundaries between consecutive components — exactly the shape `MRFDensityClusterConfig.regime_thresholds` accepts. The pre-empirical "override matches GMM auto-mode within ±2 pp at 2× ratio" sub-criterion **failed** (88.25% override vs 93.72% auto, Δ = −5.47 pp), but at low contrast (1.5×) override **wins** by +3.80 pp (75.90% vs 72.10%). The asymmetry is itself the finding (KB V14): hard `1e6` unary penalty in `_unary_from_thresholds` forfeits the Potts-smoothness benefit at high contrast (interior emitters can't be neighbor-rescued), but at low contrast it prevents the Potts amplification of GMM degeneracy that auto-mode falls into (Round 013 V13 mechanism). **Round 013 (density-ratio sweep)**: kNN-MRF accuracy 35.4% → 96.2% across {1.2, 1.5, 2.0, 3.0, 5.0}× ratios; operational floor ratio ≥ 1.65× for 75% gate, ≥ 1.85× for 85% (KB V13). Below 2×: voronoi-GMM is the better backend; with calibration data, **regime_thresholds override** is also viable for the low-contrast regime per V14. Round 012 kNN density estimator (`density_estimator=:knn`, `density_k=20`) remains the post-V11 default for high-contrast workflows; KB V12 records the patch-interior FN-band collapse from 22.5% → 6.5%. Test suite: fast tier 172/172 in ~32 s (default `Pkg.test()`), thorough tier 1002/1002 in ~36 s (gated by `SMLM_TEST_FULL=true|1|yes`) — cross-package convention shared with SMLMAnalysis / SMLMBaGoL / SMLMDriftCorrection. Remaining diagnosis priorities: patch-size sweep (P4), soft-unary override extension (P9 — new). Round dispatch is on the fork-based architecture (no tmux launcher, no lock file); orchestrator runs `/loop 60m /dispatch-round` and posts a digest + Slack notification (channel `$SLACK_CHANNEL_ID_LL_CC`) after each round. Pending external action: Priority 5 — Keith creating `JuliaSMLM/SMLMClustering.jl` so @analysis can `[sources]`-wire.

---

## Active Threads

<!-- Work explicitly in flight across multiple rounds. Each thread is a line: name, current state, owning round. Keep short. -->

(none — Task A + Task B both landed 2026-04-27)

---

## Future Priorities

<!-- Ordered list. Severity tags [CRITICAL] | [HIGH] | [MEDIUM] | [LOW] optional. Status tags TODO | IN PROGRESS | BLOCKED | DONE. -->
<!-- Human seeds these on init. Rounds reorder / update status / add discovered items, never silently drop. -->

1. [HIGH] Diagnose MRF interior false-negative band on the synthetic A431-mimic — DONE (Round 012). kNN density estimator added as `MRFDensityClusterConfig.density_estimator=:knn` mitigation. Voronoi 79.85% → kNN 89.11% headline; interior FN-rate 22.49% → 6.45%. KB V12. Operational regime bound discovered during test tuning: kNN ball must fit inside structure half-width; carry forward as input to Priority 3 (patch-size sweep).
2. [HIGH] Density-ratio operational-floor sweep — DONE (Round 013). kNN-MRF: 35.4% / 69.0% / 89.1% / 95.2% / 96.2% across ratios {1.2, 1.5, 2.0, 3.0, 5.0}×; voronoi-GMM crossover near 2× ratio. Operational floor: 75% gate at ratio ≥ 1.65×, 85% gate at ratio ≥ 1.85×, clean dominance at ratio ≥ 2.0×. Below 2×: voronoi-GMM is the better backend (Potts prior overpowers weak GMM signal). KB V13. Sweep script + 4 categorical TP/TN/FP/FN circle-plot grids in `dev/scripts/output/`.
3. [HIGH] regime_thresholds override pre-conditioning workflow — DONE (Round 014). `calibrate_regime_thresholds(smld; n_regimes=2, density_estimator=:knn, density_k=20)` helper exported from `src/backends/mrf_density.jl`; computes per-emitter log densities (Voronoi or kNN), fits 1D GMM via `_gmm_em_1d`, returns analytic Bayes decision boundaries between consecutive (mean-sorted) components — `Vector{Float64}` of length `n_regimes − 1`. ±2 pp gate at 2× ratio FAILED: 88.25% override vs 93.72% auto (Δ = −5.47 pp, hard binning forfeits Potts-smoothness benefit). At 1.5× ratio override WINS: 75.90% vs 72.10% (+3.80 pp, override avoids Potts-amplified GMM degeneracy at low contrast). KB V14. New P9 below for soft-unary override extension.
4. [HIGH] Patch-size scaling sweep. Synthetic uses 1-3 μm patches. Real biological structures span 0.3-10 μm (single hexabodies up through fiber bundles). MRF's effective patch-size range is unknown; if it fails at small patches the whole monomer-detection use case is at risk.
   - add: `dev/scripts/patch_size_sweep.jl` regenerates the synthetic at patch-size scales {0.3, 0.5, 1, 2, 5, 10 μm}, runs MRF, writes `dev/scripts/output/patch_size_sweep.csv` + `dev/scripts/output/patch_size_sweep.png`.
   - test: Round must Read the PNG and verify (a) curve is interpretable (expected: dip at ≤0.5 μm where patches drop below Voronoi neighborhood scale, then plateau), (b) MRF holds >75% accuracy in the 0.5-5 μm range. Per phase-3-extras.
   - doc: KB V14 records patch-size operating range. README MRF "When MRF works" subsection extended. — TODO
5. [MEDIUM] [BLOCKED] Push SMLMClustering to the `JuliaSMLM` GitHub org so SMLMAnalysis can pull it in via `[sources]` with `rev="main"`. @analysis's Q2 answer confirms the integration mode; @analysis is waiting on the URL + branch name. Blocked on Keith creating the repo under the org (external action). — BLOCKED
6. [LOW] MRF v2 enhancements (held for after the diagnosis loop closes): (a) graph-cuts MAP inference (binary or α-expansion) replacing ICM — `inference=:graph_cuts` symbol slot reserved but raises `ArgumentError`; needs new dep (BoykovKolmogorov.jl, GraphsFlows.jl) or in-house max-flow. (b) Soft-posterior output mode — `extras[:posterior_per_emitter]::Matrix{Float64}` (n × n_regimes) when a flag is set. (c) k-NN density estimator — landed Round 012 (`density_estimator=:knn`); this leg is closed. All remaining additive — current pipeline stays default.
   - add: TBD per (a)/(b) selection
   - test: TBD
   - doc: TBD — TODO
7. [LOW] Ratio-aware backend selector OR README guidance for low-contrast datasets. Round 013 finding: voronoi-GMM beats kNN-MRF at density ratio < 2×; MRF's Potts smoothness amplifies a weak GMM signal into uniform misclassification. Two paths: (a) document the ratio-band recommendation in README and let users pick, (b) add a meta-config that runs both and picks by some auto-criterion. Path (a) is much cheaper.
   - add: README MRF subsection extension OR `RatioAwareDensityClusterConfig` wrapper
   - test: TBD
   - doc: TBD — TODO
8. [LOW] kNN k-sensitivity sweep at fixed ratio. Round 012 surfaced the kNN-ball-radius vs structure-half-width bound; Round 013 fixed `density_k=20`. Need a sweep over k ∈ {5, 10, 20, 40, 80} at the 2× synthetic to map the (k, structure-size) operating surface. Pairs naturally with Priority 4 (patch-size sweep).
   - add: `dev/scripts/k_sensitivity_sweep.jl` — single-ratio (2×), sweep k, per-emitter accuracy + visual.
   - test: PNG read; expected curve shape: peak at moderate k, degraded at very low or very high k.
   - doc: KB entry recording optimal k for the standard synthetic, and the rule-of-thumb relating k to local density. — TODO
9. [HIGH] Soft-unary override mode for `MRFDensityClusterConfig`. Round 014 finding: the existing `regime_thresholds` override hard-pins unary cost to `0`/`1e6`, dwarfing the Potts pairwise term and forfeiting the smoothness benefit that GMM auto-mode relies on at high contrast. A soft-unary override that lets users pass calibrated `(means, vars, weights)` for each regime would build the same kind of soft Gaussian unaries (`U[i, k] = -log(w_k * N(log_rho[i] | μ_k, σ_k²))`) that GMM auto-mode produces — preserving the smoothness benefit at high contrast AND avoiding the per-ROI EM degeneracy at low contrast. Strict-Pareto improvement over both current paths.
   - add: `MRFDensityClusterConfig.regime_gaussians::Union{Nothing, NamedTuple{(:means, :vars, :weights), ...}} = nothing` field; when set, bypasses both GMM EM and hard-binning thresholds, builds soft-Gaussian unaries directly. Companion helper `calibrate_regime_gaussians(smld; ...)` returning the NamedTuple. Bonus: `calibrate_regime_thresholds` becomes a thin wrapper around `calibrate_regime_gaussians` followed by `_gaussian_decision_boundary` calls.
   - test: Thorough-tier test that the soft-unary override on the 2× synthetic matches GMM auto-mode within ±2 pp accuracy AND beats hard-threshold override at low contrast (1.5×). Companion sweep script `dev/scripts/regime_gaussians_calibration.jl` mirroring the Round 014 calibration script with the new helper.
   - doc: MRF docstring + README subsection update. KB entry follows the empirical sweep result. — TODO

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
| 012 | MRF kNN density mitigation | opus | done | kNN density estimator (`density_estimator=:knn`, `density_k=20`) added as alternative to Voronoi. Synthetic A431-mimic: 79.85% → 89.11% headline (+9.26 pp), interior FN-rate 22.49% → 6.45% (-16.04 pp). Visual: FN-band concentrated at patch interiors collapsed to boundary. Operational regime bound discovered: kNN ball radius (~√(k/πρ)) must be smaller than structure half-width or GMM split flips. Fast tier 167/167; thorough 989/989. KB V12. |
| 013 | Density-ratio operational-floor sweep | opus | done | kNN-MRF accuracy across {1.2, 1.5, 2.0, 3.0, 5.0}× ratio: 35.4 / 69.0 / 89.1 / 95.2 / 96.2%. Operational floor: 75% gate at ratio ≥ 1.65×, 85% gate at ratio ≥ 1.85×, MRF dominance at ratio ≥ 2.0×. Surprise: voronoi-GMM beats kNN-MRF below ratio 2× (Potts smoothness overpowers weak GMM signal). KB V13. New follow-up priorities P7 (low-contrast guidance) + P8 (k-sensitivity sweep). |
| 014 | regime_thresholds calibration workflow + helper API | opus | done | `calibrate_regime_thresholds(smld; n_regimes, density_estimator, density_k)` helper exported; computes Bayes decision boundaries between consecutive GMM components on per-emitter log density. ±2 pp gate at 2× ratio FAILED (88.25% override vs 93.72% auto, Δ = −5.47 pp — hard binning forfeits Potts-smoothness benefit) but at 1.5× override WINS (75.90% vs 72.10%, +3.80 pp — avoids Potts-amplified GMM degeneracy). Asymmetric outcome is the finding. Fast 172/172, thorough 1002/1002. KB V14. New P9 for soft-unary override extension (strict-Pareto fix). |

---

## Stop Conditions

<!-- Halt Phase 2 before any work if any of these are met. Add a line to halt autonomous operation without touching the lock. -->

- (none — diagnosis loop active; Priority 5 is BLOCKED on external action but does not gate other rounds.)

<!-- Examples:
- `Test suite broken on main branch (block all rounds until fixed)`
- `Paused by kalidke on 2026-04-15: waiting for external review`
- `All priorities DONE — project complete`
-->
