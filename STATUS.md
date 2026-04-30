# STATUS — SMLMClustering

*Live project state. Read in Phase 1 of every round. Updated in Phase 4.*

---

## Current State

All five labeling backends (`DBSCANConfig`, `HDBSCANConfig`, `HierarchicalConfig`, `VoronoiConfig`, `MRFDensityClusterConfig`) and both spatial-statistic backends (`HopkinsConfig`, `VoronoiDensityConfig`) are live end-to-end with full documentation, README, and api_overview. `cluster(smld, cfg)` is non-mutating (KB V9 deepcopy); `cluster_statistics(smld, cfg)` is pass-through (KB V10 same-reference). **Round 016 (kNN k-sensitivity sweep) closed (KB V17)**: at fixed 2× density ratio + `patch_scale=0.75` (~1.5 μm patches), the kNN-MRF accuracy curve over `density_k ∈ {5, 10, 15, 20, 30, 50, 80}` is concave-up with peak at k=15 (91.23%) and k=20 within rounding (91.01%); all 7 k-values clear the 75% gate, six (k=10..80) clear the 85% gate. Precision/recall asymmetry is the underlying mechanism: low k → low precision via FP-dominant smoothing-amplified noise (k=5 prec=78.83%, FP=1295), high k → low recall via FN-dominant density-mean regression (k=80 rec=70.35%, FN=1652). The V12-default k=20 is reaffirmed; the marginal k=15 advantage (+0.22 pp) is within rounding. With V12 / V13 / V16 / V17 the (ratio, patch-size, k) operating surface is mapped at three orthogonal slices through the 2×-ratio-1.0/1.5-μm-patches-k=20 anchor. **Round 015 (patch-size scaling sweep) closed (KB V16)**: kNN-MRF at 2× contrast holds 88-89% across 1.0-3.0 μm patches, climbs through the 75% gate at 0.5 μm and the 85% gate between 0.5 and 1.0 μm; collapses to 62.7% at 0.3 μm (V12 regime bound). voronoi-GMM is a strict-loss baseline at 2× ratio. **Soft-unary calibration landed direct 2026-04-30 (KB V15)**: `MRFDensityClusterConfig.regime_gaussians` + `calibrate_regime_gaussians(...)` bypass per-ROI EM while preserving soft Gaussian unaries; held-out 2× soft 92.66% vs auto 93.72% / hard 88.25%; 1.5× soft 76.95% vs hard 75.90% / auto 72.10%. **Round 014 (regime_thresholds)**: hard binning loses ~5 pp at high contrast (forfeits Potts neighbor rescue) but wins ~4 pp at 1.5× contrast (avoids Potts-amplified GMM degeneracy) — KB V14. **Round 013 (density-ratio sweep)**: kNN-MRF 35.4 → 96.2% across {1.2, 1.5, 2.0, 3.0, 5.0}× — KB V13, operational floor ratio ≥ 1.65× / ≥ 1.85× / ≥ 2.0× for the 75 / 85 / clean-dominance bars. Round 012 kNN-density default (k=20) — KB V12 patch-interior FN-band collapse 22.5 → 6.5%. Test suite: fast tier 185/185 in ~33 s, thorough tier 1025/1025 in ~38 s (gated by `SMLM_TEST_FULL=true|1|yes`). Round dispatch is on the fork-based architecture (no tmux launcher, no lock file); orchestrator runs `/loop 60m /dispatch-round` and posts a digest + Slack notification (channel `$SLACK_CHANNEL_ID_LL_CC`) after each round. Pending external action: Priority 5 — Keith creating `JuliaSMLM/SMLMClustering.jl` so @analysis can `[sources]`-wire. P6/P7 remain LOW-priority with TBD trio specs (Q6/Q7 OPEN).

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
3. [HIGH] regime_thresholds override pre-conditioning workflow — DONE (Round 014). `calibrate_regime_thresholds(smld; n_regimes=2, density_estimator=:knn, density_k=20)` helper exported from `src/backends/mrf_density.jl`; computes per-emitter log densities (Voronoi or kNN), fits 1D GMM via `_gmm_em_1d`, returns analytic Bayes decision boundaries between consecutive (mean-sorted) components — `Vector{Float64}` of length `n_regimes − 1`. ±2 pp gate at 2× ratio FAILED: 88.25% override vs 93.72% auto (Δ = −5.47 pp, hard binning forfeits Potts-smoothness benefit). At 1.5× ratio override WINS: 75.90% vs 72.10% (+3.80 pp, override avoids Potts-amplified GMM degeneracy at low contrast). KB V14. Soft-unary follow-up landed as Priority 9 / KB V15.
4. [HIGH] Patch-size scaling sweep — DONE (Round 015). `dev/scripts/patch_size_sweep.jl` sweeps `patch_scale ∈ {0.15, 0.25, 0.5, 1.0, 1.5}` (nominal mean rect length `{0.3, 0.5, 1.0, 2.0, 3.0} μm`) at fixed 2× density ratio with kNN-MRF (k=20) + voronoi-GMM baseline. kNN-MRF accuracy: 62.67% / 79.25% / 88.72% / 89.11% / 88.31% — passes 75% gate at ≥0.5 μm, 85% gate at ≥1.0 μm, 88-89% plateau across 1.0-3.0 μm. Failure at 0.3 μm matches V12's predicted regime bound (kNN ball radius ~80 nm comparable to structure half-width 75-150 nm — ball spills into background). Sweep truncated at 3.0 μm because 5×5 μm ROI cannot fit larger patches with non-overlap; KB V16. Visual diagnostics: `dev/scripts/output/patch_size_sweep.png` (accuracy curve, log-x) + `dev/scripts/output/round_015_smlmrender_categorical.png` (5-panel SMLMRender categorical TP/TN/FP/FN render, kNN-MRF predictions only).
5. [MEDIUM] [BLOCKED] Push SMLMClustering to the `JuliaSMLM` GitHub org so SMLMAnalysis can pull it in via `[sources]` with `rev="main"`. @analysis's Q2 answer confirms the integration mode; @analysis is waiting on the URL + branch name. Blocked on Keith creating the repo under the org (external action). — BLOCKED
6. [LOW] MRF v2 enhancements (held for after the diagnosis loop closes): (a) graph-cuts MAP inference (binary or α-expansion) replacing ICM — `inference=:graph_cuts` symbol slot reserved but raises `ArgumentError`; needs new dep (BoykovKolmogorov.jl, GraphsFlows.jl) or in-house max-flow. (b) Soft-posterior output mode — `extras[:posterior_per_emitter]::Matrix{Float64}` (n × n_regimes) when a flag is set. (c) k-NN density estimator — landed Round 012 (`density_estimator=:knn`); this leg is closed. All remaining additive — current pipeline stays default.
   - add: TBD per (a)/(b) selection
   - test: TBD
   - doc: TBD — TODO
7. [LOW] Ratio-aware backend selector OR README guidance for low-contrast datasets. Round 013 finding: voronoi-GMM beats kNN-MRF at density ratio < 2×; MRF's Potts smoothness amplifies a weak GMM signal into uniform misclassification. Two paths: (a) document the ratio-band recommendation in README and let users pick, (b) add a meta-config that runs both and picks by some auto-criterion. Path (a) is much cheaper.
   - add: README MRF subsection extension OR `RatioAwareDensityClusterConfig` wrapper
   - test: TBD
   - doc: TBD — TODO
8. [LOW] kNN k-sensitivity sweep at fixed ratio — DONE (Round 016). `dev/scripts/k_sensitivity_sweep.jl` sweeps `density_k ∈ {5, 10, 15, 20, 30, 50, 80}` at fixed 2× density ratio + `patch_scale=0.75` (~1.5 μm patches). Concave-up curve: 85.15 / 89.31 / 91.23 / 91.01 / 90.50 / 88.26 / 87.32%. Peak at k=15-20 (precision/recall crossover), all 7 values clear 75% gate, six (k=10..80) clear 85% gate. V12-default k=20 reaffirmed (marginal k=15 advantage +0.22 pp within rounding). KB V17 records the bound: r_k(k=20)=80 nm at ρ_high=1000/μm² is the right ball-radius scale for ~1.5 μm patches with 75-150 nm rect-thin half-widths. Visual: `dev/scripts/output/round_016_smlmrender_categorical.png` — 7-panel SMLMRender categorical TP/TN/FP/FN, FP-green at low k → FN-red at high k. Sweep at single (ratio, patch-size); (k, ratio) and (k, patch-size) joint surfaces remain open.
9. [HIGH] Soft-unary override mode for `MRFDensityClusterConfig` — DONE direct 2026-04-30. `MRFDensityClusterConfig.regime_gaussians::Union{Nothing, NamedTuple{(:means, :vars, :weights)}} = nothing` field added; when set, bypasses both GMM EM and hard-binning thresholds and builds soft Gaussian unaries directly. `calibrate_regime_gaussians(smld; ...)` returns the NamedTuple; `calibrate_regime_thresholds` is now a thin wrapper around that helper plus `_gaussian_decision_boundary`. Fast regression proves soft unaries allow neighbor rescue where hard thresholds cannot; thorough tests cover helper shape and end-to-end use. Held-out A431-mimic check: 2× soft 92.66% vs auto 93.72% (within ±2 pp) and hard 88.25%; 1.5× soft 76.95% vs hard 75.90% and auto 72.10%. README/api_overview/docstrings updated. KB V15.

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
| 015 | Patch-size scaling sweep | sonnet | done | Patch-size sweep at fixed 2× ratio, kNN-MRF (k=20) + voronoi-GMM baseline across nominal patch sizes `{0.3, 0.5, 1.0, 2.0, 3.0} μm`. kNN-MRF: 62.7% / 79.3% / 88.7% / 89.1% / 88.3% — 75% gate at ≥0.5 μm, 85% gate at ≥1.0 μm, plateau 88-89% across 1.0-3.0 μm; collapses at 0.3 μm (kNN ball spills into background, V12 regime bound). voronoi-GMM strict-loss across this contrast (~65-72%). 1.0 μm reproduces Round 012's 89.11% headline (positive control). Sweep ROI-truncated at 3.0 μm. KB V16. Bundles uncommitted V15 source/test/doc changes (regime_gaussians soft-unary). |
| 016 | kNN k-sensitivity sweep | sonnet | done | k-sweep at fixed 2× ratio + ~1.5 μm patches, kNN-MRF across k ∈ {5, 10, 15, 20, 30, 50, 80}: 85.15 / 89.31 / 91.23 / 91.01 / 90.50 / 88.26 / 87.32%. Concave-up curve, peak at k=15-20 (precision/recall crossover); all clear 75% gate, six clear 85%. V12-default k=20 reaffirmed (Δ vs k=15 peak = +0.22 pp, within rounding). Precision/recall asymmetry: low k → FP-dominant (smoothing-amplified noise), high k → FN-dominant (density-mean regression). KB V17. Two OPEN questions Q6/Q7 added for P6/P7 trio specs. |

---

## Stop Conditions

<!-- Halt Phase 2 before any work if any of these are met. Add a line to halt autonomous operation without touching the lock. -->

- (none — diagnosis loop active; Priority 5 is BLOCKED on external action but does not gate other rounds.)

<!-- Examples:
- `Test suite broken on main branch (block all rounds until fixed)`
- `Paused by kalidke on 2026-04-15: waiting for external review`
- `All priorities DONE — project complete`
-->
