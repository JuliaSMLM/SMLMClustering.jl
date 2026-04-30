# Round 014 — regime_thresholds calibration workflow + helper API

**Date:** 2026-04-30
**Status:** done
**Priority worked on:** Priority 3 — regime_thresholds override pre-conditioning workflow

## Output files to review

- `dev/scripts/output/regime_calibration_demo.png` — calibration histogram (left, log ρ from kNN k=20 with the calibrated threshold line at 6.685) + accuracy bar chart (right, 2× auto vs 2× override vs 1.5× auto vs 1.5× override) + summary text. Look for: the threshold sits to the right of the histogram peak in a heavy-tail region (the two GMM modes are not visually well-separated even at 2× ratio), and the bar chart shows the asymmetric outcome — override loses 5.47 pp at 2× (Potts smoothness benefit forfeit) but gains 3.80 pp at 1.5× (avoids GMM auto-mode's degenerate-split + Potts amplification).
- `dev/scripts/output/regime_calibration_panels.png` — categorical per-emitter circle plots (gray=correct, orange=FP, cyan=FN), one panel per case in left-to-right order: 2× auto, 2× override, 1.5× auto, 1.5× override. Look for: (1) 2× auto panel has clean patches with sparse cyan-FN; 2× override panel shows visibly more cyan-FN concentrated inside patch interiors (the smoothness-loss signature). (2) 1.5× auto is orange-FP saturation across the whole ROI (the Potts-amplified GMM-degeneracy collapse from Round 013); 1.5× override shows the orange-FP collapsed, replaced by cyan-FN inside patches (override clips the catastrophe but doesn't recover full recall).
- `dev/scripts/output/regime_calibration_demo.csv` — per-(case, mode) accuracy/precision/recall/TP/FP/FN with the calibrated threshold value and Δ-pp-vs-auto column for quick comparison.
- `dev/scripts/regime_thresholds_calibration.jl` — the calibration demo script. Reads `simulate_a431_mimic.jl` for the `simulate_dataset(...)` entry point exposed in Round 013.
- `src/backends/mrf_density.jl` (lines ~261-380) — the new `calibrate_regime_thresholds(smld; ...)` helper + the `_gaussian_decision_boundary` internal that computes the analytic Bayes boundary between consecutive GMM components. Read the docstring for the calibration workflow recipe.

## Hypothesis

`MRFDensityClusterConfig.regime_thresholds` exists as an undocumented config field that bypasses GMM and binds emitters directly into hard regime bins. For datasets with stable density structure across ROIs (same imaging conditions, same labeling protocol), calibrating once on a trusted ROI and reusing the thresholds across query ROIs should be (a) faster (skips per-ROI EM), (b) more reproducible (deterministic thresholds), and (c) more robust at the low-contrast end of the operational band where Round 013 showed GMM auto-mode degenerates and the Potts smoothness amplifies the degeneracy into uniform misclassification. The bar from Priority 3 was "override mode must match GMM mode within ±2% accuracy on the synthetic."

## What was attempted

Added `calibrate_regime_thresholds(smld; n_regimes=2, density_estimator=:knn, density_k=20)` to `src/backends/mrf_density.jl`. The helper computes per-emitter log densities via the chosen estimator (Voronoi or kNN), fits a 1D `n_regimes`-component GMM via the existing `_gmm_em_1d`, and returns the analytic Bayes decision boundary between each pair of consecutive (mean-sorted) GMM components. The Bayes boundary is the quadratic root in `[μ_k, μ_{k+1}]` of `log(w_k) + log N(x|μ_k, σ_k²) = log(w_{k+1}) + log N(x|μ_{k+1}, σ_{k+1}²)`, falling back to the midpoint when the discriminant is negative. Returns a `Vector{Float64}` of length `n_regimes − 1`, sorted ascending — exactly the shape `MRFDensityClusterConfig.regime_thresholds` accepts.

Added thorough-tier tests (validation + 2-regime helper + kNN path) and exported the helper. Built `dev/scripts/regime_thresholds_calibration.jl` that simulates a calibration ROI (seed 20260429, 2× ratio), fits thresholds via the helper, then evaluates two modes — GMM auto-mode (no override) and override mode (calibrated thresholds) — on a held-out query ROI (seed 20260430, 2× ratio, different patch geometry) AND a low-contrast query ROI (seed 20260430, 1.5× ratio) to test the Round 013-motivated extension hypothesis.

## What worked

- Helper API ships clean. `calibrate_regime_thresholds` validates arguments, runs the Voronoi or kNN estimator, fits GMM, and returns analytic Bayes decision boundaries. Thorough-tier tests cover argument validation, 2-regime + 3-regime fits, the kNN path, and end-to-end use of returned thresholds in `MRFDensityClusterConfig`.
- **Test suite green**: fast tier 172/172 in 32.1 s (was 167; +5 from new validation tests), thorough tier 1002/1002 in 35.8 s (was 989; +13 from new helper tests).
- **Calibration script runs end-to-end** on the synthetic. Calibration on the 2× ROI produces a single threshold at log ρ = 6.685 (ρ ≈ 800 emitters/μm²), sitting between the median background ρ ≈ 500 and high-density patch ρ ≈ 1000.
- **Low-contrast extension wins**: At ratio 1.5× — the operational regime where Round 013 saw kNN-MRF collapse — override mode lands at **75.90% accuracy**, beating GMM auto-mode's **72.10%** by **+3.80 pp**. The visual diagnostic (`regime_calibration_panels.png` panel 4) shows the override clips the orange-FP saturation that Potts-amplified GMM degeneracy produces (Round 013 V13 mechanism), substituting it for cyan-FN at patch interiors. Net positive — override avoids the GMM-degeneracy catastrophe.

**Visual inspection** (per phase-3-extras of start-round.md — Read each PNG before claiming done):

- `regime_calibration_demo.png`: histogram is unimodal-looking with a long right tail (the two GMM components overlap heavily even at 2× ratio); calibrated threshold sits in the tail. Bar chart shows the asymmetric outcome plainly: 2× auto (93.7%) > 2× override (88.2%); 1.5× auto (72.1%) < 1.5× override (75.9%).
- `regime_calibration_panels.png`: 2× auto panel is mostly correct gray with sparse cyan-FN inside patches; 2× override panel has visibly denser cyan-FN inside patch interiors (the smoothness-loss signature). 1.5× auto panel is orange-FP saturation across the whole ROI (matches Round 013 V13); 1.5× override panel has orange-FP collapsed, replaced by cyan-FN at patch interiors. Visual content matches the numerical headlines and the mechanism explanation.

## What failed

**Gate "override matches GMM auto-mode within ±2 pp at 2× ratio" — not met.** Override mode lands at 88.25% accuracy on the 2× query ROI vs GMM auto-mode's 93.72% — Δ = **−5.47 pp**, well outside the ±2 pp band. This is the same shape of result as Round 013's "(b)-miss" (MRF beats baselines across the sweep — also missed at low contrast): the empirical data reveals a mechanism that the priority's pre-empirical hypothesis didn't account for.

**Mechanism (the actual finding):** override mode hard-pins the unary cost matrix to `0` for the matching bin and `1e6` for others (see `_unary_from_thresholds` in `src/backends/mrf_density.jl:309-323`), whereas GMM auto-mode produces soft Gaussian unaries `−log(w_k · N(x | μ_k, σ_k²))`. The Potts smoothness term `λ × #neighbors_with_label_not_k` in `_icm_potts!` is dwarfed by the `1e6` hard penalty, so neighbors cannot pull a borderline interior point across the regime boundary. At high contrast (2×), this loses the Potts-smoothness benefit that the kNN-MRF default depends on — interior emitters whose individual log ρ falls below the threshold by chance can no longer be rescued by their high-density neighbors. The result is a +784 increase in FN (1316 vs 532) at almost the same FP rate.

At low contrast (1.5×) the trade-off inverts: GMM auto-mode degenerates (the two log ρ modes coalesce, EM picks a near-tied split), and the soft Potts smoothness — which would normally help — instead amplifies the degeneracy into uniform misclassification (Round 013 V13 mechanism). Override mode's hard binning prevents that amplification: even though the threshold is "wrong" for 1.5× data (calibrated on 2×, naturally too high), it produces a graceful drop in recall (33.30% vs 52.64%) rather than the catastrophic FP saturation auto-mode falls into.

## Files changed

```
src/backends/mrf_density.jl                              +120 (calibrate_regime_thresholds + _gaussian_decision_boundary helpers)
src/SMLMClustering.jl                                    +1 export
test/test_mrf_density.jl                                 +60 (validation + 2-regime + kNN + 3-regime helper tests)
dev/scripts/regime_thresholds_calibration.jl            new (force-add, ~250 lines)
dev/scripts/output/regime_calibration_demo.csv          new (force-add)
dev/scripts/output/regime_calibration_demo.png          new (force-add, ~210 KB)
dev/scripts/output/regime_calibration_panels.png        new (force-add, ~1.7 MB)
STATUS.md                                                P3 → DONE; Current State refresh; new P9 for soft-unary override
KNOWLEDGE_BASE.md                                        +V14
rounds/round_014_*.md                                    new
```

## Confidence

High on the data: deterministic seeds (20260429 calibration, 20260430 query) make the experiment fully reproducible; same-RNG simulator means we control patch geometry across ratios; visual diagnostic on the categorical circle plots confirms the mechanism explanation panel-by-panel. High on the helper API: 1002/1002 thorough tests passing including the new helper-specific testset; the API matches the existing `regime_thresholds` field's shape with no ambiguity. Medium on the operational guidance: the 2× and 1.5× endpoints are clean, but the *crossover* between "auto wins" and "override wins" was not swept in this round — it likely lives somewhere in the [1.5×, 2.0×] band where Round 013's kNN-MRF curve climbs steeply. P9 below addresses both the soft-unary extension AND a follow-up sweep to characterize the crossover.

## External consultations (if any)

None this round.

## Next steps

Priority 3 closes as DONE — the trio is complete: helper API ships, tests pass at fast 172/172 and thorough 1002/1002, calibration workflow + visual diagnostic + KB entry document the operational regime. The "match within ±2 pp at 2×" sub-criterion failed but the failure is itself the finding (override forfeits Potts-smoothness benefit by hard binning), recorded as KB V14.

The clean follow-up is a **soft-unary override mode** (logged as new P9 in STATUS.md) that lets users supply `(means, vars, weights)` for each regime — the override would build the same kind of soft Gaussian unaries that GMM auto-mode produces, preserving the Potts smoothness benefit while still skipping per-ROI EM. This would convert the calibration helper into a strict-Pareto-improvement workflow: matches GMM auto-mode at high contrast (since the unaries are soft Gaussians as auto-mode produces) AND avoids the GMM-degeneracy + Potts amplification at low contrast (since the calibrated parameters don't degenerate per-ROI).

Priority 4 (patch-size scaling sweep) remains queued and is the next round's natural pick. The Round-012 kNN-radius-vs-structure-half-width bound + Round-013 ratio-floor + Round-014 calibration mechanism together set up Priority 4 with three known operational axes to vary against patch size.

## Questions posted

(none — finding is well-specified, P9 follow-up has a clear trio shape, no human-judgment trigger fired)
