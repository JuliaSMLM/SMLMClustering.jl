# Round 013 — Density-ratio operational-floor sweep

**Date:** 2026-04-30
**Status:** done
**Priority worked on:** Priority 2 — Density-ratio operational-floor sweep

## Output files to review

- `dev/scripts/output/density_ratio_sweep.png` — accuracy-vs-ratio curve, all 4 backends, with 75%/85% gate lines. Look for: kNN-MRF crossover with voronoi-GMM near 2× ratio.
- `dev/scripts/output/density_ratio_sweep.csv` — per-(backend, ratio) accuracy/precision/recall/TP/FP/FN.
- `dev/scripts/output/density_ratio_sweep_mrf_knn_panels.png` — categorical TP/TN/FP/FN per-emitter circle plot, 5 ratio panels for kNN-MRF. Look for: massive orange-FP saturation at 1.2× collapsing to clean classification at ≥3×.
- `dev/scripts/output/density_ratio_sweep_voronoi_gmm_panels.png` — same shape for voronoi-GMM. Look for: more graceful degradation at low ratio (it's the low-ratio winner).
- `dev/scripts/output/density_ratio_sweep_dbscan_panels.png` — same for DBSCAN. Look for: 100% recall at the cost of catastrophic FP — DBSCAN labels nearly everything "high".
- `dev/scripts/output/density_ratio_sweep_hdbscan_panels.png` — same for HDBSCAN. Look for: flat ~55-65% accuracy — its parameter regime doesn't track ratio.
- `dev/scripts/density_ratio_sweep.jl` — the sweep script. Reads sim helpers via `include("simulate_a431_mimic.jl")`.
- `dev/scripts/simulate_a431_mimic.jl` — refactored to expose `simulate_dataset(; rho_low, rho_high_bonus, seed, verbose)` returning a NamedTuple. `main()` is gated behind `abspath(PROGRAM_FILE) == @__FILE__` so `include()` is side-effect-free.

## Hypothesis

A 5×5 μm A431-mimic regenerated at controlled density ratios — keeping patch geometry fixed via the same `SIM_SEED`, varying `RHO_HIGH_BONUS` so high-density = `RHO_LOW × ratio` — should let us trace MRF's accuracy floor as the density contrast shrinks. Expected qualitative shape: monotonic rise with ratio for kNN-MRF, dominance over baselines across the sweep, and a clean operational floor (min ratio at which kNN-MRF accuracy ≥ 75%) somewhere between 1.2× and 2×. Real dSTORM data ratios run 1.2–3× depending on labeling efficiency and probe density, so the sweep informs deployment decisions for the genmab paper push.

## What was attempted

Refactored `dev/scripts/simulate_a431_mimic.jl` to expose `simulate_dataset(; rho_low, rho_high_bonus, seed, verbose)` (returning `(smld, ground_truth, patches, stats)` NamedTuple) without rerunning the default save-to-jld2 main path on `include()`. Same `SIM_SEED` produces identical patch geometry across calls — the experiment is controlled, only the per-patch emitter density varies.

Built `dev/scripts/density_ratio_sweep.jl` that sweeps ratios {1.2×, 1.5×, 2.0×, 3.0×, 5.0×}. For each ratio it runs four backends and records per-emitter (high vs low) accuracy/precision/recall:

1. **kNN-MRF** (`MRFDensityClusterConfig(n_regimes=2, density_estimator=:knn, density_k=20)`) — post-Round-012 winner.
2. **DBSCAN** (`DBSCANConfig(eps_nm=100.0, min_points=5)`) — noise → low, cluster → high.
3. **HDBSCAN** (`HDBSCANConfig(min_points=5, knn_graph_k=30)`) — same noise/cluster mapping.
4. **Voronoi-GMM** — `VoronoiDensityConfig` density vector → 2-component 1D GMM EM (no MRF smoothing) — baseline for "what does the GMM fit alone tell you?".

Outputs: a CSV with per-(backend, ratio) rows; a 4-curve aggregate accuracy plot with the 75% and 85% gate lines; and per-backend categorical TP/TN/FP/FN circle-plot grids (one figure per backend, 5 ratio panels per figure). All in `dev/scripts/output/`.

The sweep adds 9 force-tracked files to git (`git add -f`); `dev/` is gitignored at the project root.

## What worked

Headline numbers (mrf_knn = kNN k=20, the post-Round-012 default):

| Ratio  | kNN-MRF | DBSCAN | HDBSCAN | Voronoi-GMM |
|--------|---------|--------|---------|-------------|
| 1.2×   | 35.4%   | 29.3%  | 49.4%   | **64.5%**   |
| 1.5×   | 69.0%   | 34.6%  | 53.0%   | 66.5%       |
| 2.0×   | **89.1%**| 41.3% | 57.2%   | 71.8%       |
| 3.0×   | **95.2%**| 51.6% | 65.3%   | 69.1%       |
| 5.0×   | **96.2%**| 64.1% | 63.0%   | 86.8%       |

(Bold: per-row best.)

Operational floor (kNN-MRF accuracy ≥ 75%) by linear interpolation: ratio ≈ **1.65×**. The 85% gate falls at ratio ≈ **1.85×**. Sanity check: the 2.0× column reproduces Round 012's 89.11% headline exactly (same patches, same RNG, same backend) — the refactor is correct.

**Visual inspection** (per phase-3-extras of start-round.md — Read each PNG before claiming done):

- `density_ratio_sweep.png`: the kNN-MRF (red) curve climbs steeply from 1.2× → 2.0× then plateaus at ~95-96% above 3×. Voronoi-GMM (yellow) is the surprise lead below ratio 2× — at 1.2× it's at 64.5% while MRF is at 35.4%. The two curves cross between 1.5× and 2.0×, where MRF takes over for good. DBSCAN (green) and HDBSCAN (purple) trail throughout; HDBSCAN is flat ~50-65% indicating its parameter regime doesn't engage with this density-segmentation task.
- `density_ratio_sweep_mrf_knn_panels.png`: at 1.2× the entire ROI is orange (FP) saturated — kNN-MRF over-calls "high" almost everywhere because the GMM regime split degenerates when the two density modes are too close. At 1.5× orange is reduced but still dominates the background; cyan FN appear inside true patches. At 2.0× the panel cleans up dramatically — light orange edges, sparse cyan inside patches, clear patch shapes. By 3.0× and 5.0× nearly all points are gray (correct); patches are saturated correctly. The visual matches the numerical headline.
- `density_ratio_sweep_voronoi_gmm_panels.png` (spot-checked): degrades more gracefully than MRF at low ratio — the 1.2× and 1.5× panels still show patch shapes through the noise rather than going entirely orange. This is the structural reason voronoi-GMM wins below 2× — without the Potts smoothness term, low-ratio noise can't get amplified by neighbor-coherence forcing.

## What failed

**Gate (b) — "MRF beats baselines across the sweep" — was not met.** kNN-MRF is the *worst* of the four backends at ratio 1.2× (35.4% vs voronoi-GMM's 64.5%) and roughly tied with voronoi-GMM at 1.5× (69.0% vs 66.5%). The MRF Potts smoothness term that helps interior coherence at high contrast becomes a liability at low contrast — it amplifies the GMM's degenerate regime split into uniform misclassification. This is a real result, not a bug: the Potts prior assumes the labels are smoothly distributed, which is only useful when the GMM has actually identified two regimes; at very low contrast the GMM means coalesce and the prior smooths everything to one label.

This is consistent with the operational regime bound discovered in Round 012 (kNN ball radius vs structure half-width). Both failure modes share a root cause: when the GMM signal is weaker than the MRF prior assumes, the prior swamps the data.

## Files changed

```
dev/scripts/density_ratio_sweep.jl                              new (force-add)
dev/scripts/simulate_a431_mimic.jl                              refactored, force-add (was untracked)
dev/scripts/output/density_ratio_sweep.csv                      new
dev/scripts/output/density_ratio_sweep.png                      new
dev/scripts/output/density_ratio_sweep_mrf_knn_panels.png       new
dev/scripts/output/density_ratio_sweep_dbscan_panels.png        new
dev/scripts/output/density_ratio_sweep_hdbscan_panels.png       new
dev/scripts/output/density_ratio_sweep_voronoi_gmm_panels.png   new
STATUS.md                                                       P2 → DONE; Current State refresh
KNOWLEDGE_BASE.md                                               +V13
README.md                                                       MRF "When it works" subsection
rounds/round_013_*.md                                           new
```

No `src/` or `test/` changes — the round is a measurement on existing code. Test suite remains 167/167 fast / 989/989 thorough as of Round 012.

## Confidence

High on the data: 5 ratios × 4 backends × ~13k-22k emitters per ratio is a fully-specified controlled experiment, the 2.0× column reproduces Round 012's headline as a sanity check, and visual inspection of the per-emitter classification panels confirms the curve. Medium on the interpretation of the (b) miss: the explanation (Potts prior overpowers a weak GMM signal at low contrast) is consistent with everything we know about MRFs and is a sensible hypothesis, but it could also be a knob-tuning issue (smaller `density_k`, smaller `smoothness_lambda`, or a non-Potts pairwise might rescue MRF at low ratio). That hypothesis goes into the next-priority backlog rather than getting tested here — it's outside the (a)/(b)/(c) scope of P2.

## External consultations (if any)

None this round.

## Next steps

Priority 2 closes as DONE. The operational floor is **ratio ≥ 1.65× for the 75% gate** and **ratio ≥ 1.85× for the 85% gate**, with clean MRF dominance over all baselines starting at **ratio ≥ 2.0×**. The (b)-miss surfaces a new follow-up: at low contrast, **voronoi-GMM is the better backend** — the round logs a new LOW-priority item to add a "ratio-aware backend selector" or document the recommendation in README.

Priority 3 (regime_thresholds calibration workflow) and Priority 4 (patch-size sweep) remain queued. Round 012's kNN-radius-vs-structure-half-width bound is directly relevant to Priority 4; both should be considered together when evaluating MRF's full operating envelope.

## Questions posted

(none — sweep produced unambiguous, inspectable results; no human-judgment trigger fired)
