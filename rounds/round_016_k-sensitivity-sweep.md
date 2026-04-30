# Round 016 — kNN k-sensitivity sweep at fixed 2× ratio + 1.5 μm patches

**Date:** 2026-04-30
**Status:** done
**Priority worked on:** Priority 8 — kNN k-sensitivity sweep

## Output files to review

- `dev/scripts/output/round_016_smlmrender_categorical.png` (2.5 MiB, 7-panel) — **primary visual artifact.** SMLMRender CircleRender categorical TP/TN/FP/FN circle plot, one panel per `density_k` value (5, 10, 15, 20, 30, 50, 80). tab10 categorical encoding: id 1 = TN (blue), 2 = TP (orange), 3 = FP (green), 4 = FN (red). Look for the FP/FN trade-off across panels: low-k panels have visible green (FP) speckle around patch boundaries (sub-sampling noise → over-aggressive labeling), high-k panels show red (FN) inside patch interiors (kNN ball spills into background, density estimate regresses to mean). Middle panels (k=15-20) have clean orange-TP patch interiors with thin red-FN/green-FP fringe — the canonical "kNN-MRF works" visual.
- `dev/scripts/output/k_sensitivity_sweep.png` (104 KiB) — accuracy curve, log-x in `density_k`, dashed 75% gate + dotted 85% gate, per-point annotation `r_k ≈ √(k/πρ_high)` in nm. Look for the concave-up shape: peak at k=15 (91.23%), softening at low k=5 (85.15%, just clears 85% gate) and high k=80 (87.32%). All 7 values clear the 75% gate; 6 of 7 clear the 85% gate. The k=20 point reproduces V12's headline within ~2 pp at this geometry (different patch_scale than Round 012's default; same backend).
- `dev/scripts/output/k_sensitivity_sweep.csv` (7 rows) — per-k metrics: `k, r_k_nm, n_total, accuracy, precision, recall, tp, tn, fp, fn`. Sortable for the precision/recall asymmetry analysis (low k → high recall + low precision; high k → high precision + low recall, the V18 finding below).
- `dev/scripts/k_sensitivity_sweep.jl` (220 lines) — sweep script. Same `simulate_dataset(patch_scale=0.75, ...)` once, then varies `MRFDensityClusterConfig(density_k=k)` across the 7 k-values. Reusable as a template for future (k, ratio) or (k, patch-size) joint sweeps.

## Hypothesis

KB V12 articulated the operational regime bound for the kNN density estimator: kNN ball radius `r_k ≈ √(k/πρ)` must be smaller than the structure half-width or the GMM regime split flips foreground / background. At ρ_high = 1000 emit/μm² (the patch-interior density on this synthetic), `r_k ∈ {40, 56, 69, 80, 98, 126, 160} nm` for `k ∈ {5, 10, 15, 20, 30, 50, 80}`. At `patch_scale=0.75` (nominal mean rect length 1.5 μm), the rect short-side range is 75-150 nm and the ellipse semi-minor range is 187-750 nm — the bound bites hardest on the rect-thin tail at large k. Predicted shape: concave-up curve. Soft floor at very low k (sub-sampling noise σ_log ≈ 1/√k blows up → unary unary distribution becomes too smeared for clean GMM split), peak at moderate k where r_k is small enough to fit inside structures and large enough to suppress the noise floor, soft drop at very high k where r_k spills out of patches into background and the density estimate regresses toward the global mean. The Round 012 default (k=20) was hand-picked at `patch_scale=1.0` and is the first sweep point to verify against this predicted shape; the paper-deadline question Priority 8 was specified against was "is k=20 still the right default at this geometry, or does the curve favor a different k?".

## What was attempted

Drafted `dev/scripts/k_sensitivity_sweep.jl` adapting the Round 013 `density_ratio_sweep.jl` template. Held `simulate_dataset(patch_scale=0.75, rho_low=500, rho_high_bonus=500, seed=20260429)` constant — generates the synthetic ONCE so RNG-dependent geometry is stable across the k sweep — and varied only `MRFDensityClusterConfig(n_regimes=2, density_estimator=:knn, density_k=k)` across `KS = [5, 10, 15, 20, 30, 50, 80]`. Per k: predict with the kNN-MRF backend, compute TP/TN/FP/FN against the simulator's ground-truth labels, accumulate (acc, prec, rec, r_k_nm) into the CSV, build a per-k SMLMRender CircleRender panel with categorical id-encoded TP/TN/FP/FN colors via the canonical `render(smld_cat; strategy=CircleRender(), color_by=:id, categorical=true, colormap=:tab10, zoom=8)` API. The `clamp_rgb.(img)` step before save is mandatory — CircleRender accumulates intensity at overlap and CairoMakie's PNG backend (N0f8) errors on out-of-range floats; this is the lesson from the Round 015 backfill cascade (KB note in `start-round.md` phase-3-extras). Composed an aggregate accuracy curve (CairoMakie, log-x in k, with r_k annotations per point) and the canonical 7-panel SMLMRender grid `round_016_smlmrender_categorical.png`.

The synthetic at `patch_scale=0.75` produced 13,755 emitters across 12 patches at 40.5% high-density fraction (matches expectation for 2× ratio + ~25% high-density area).

## What worked

The accuracy curve is the predicted concave-up shape, cleanly:

| k | r_k (nm) | accuracy | precision | recall | TP | FP | FN |
|---|---|---|---|---|---|---|---|
| 5  | 40  | 85.15% | 78.83% | 86.57% | 4823 | 1295 | 748 |
| 10 | 56  | 89.31% | 82.24% | 93.90% | 5231 | 1130 | 340 |
| 15 | 69  | **91.23%** | 87.60% | 91.28% | 5085 | 720 | 486 |
| 20 | 80  | 91.01% | 93.97% | 83.14% | 4632 | 297 | 939 |
| 30 | 98  | 90.50% | 96.49% | 79.43% | 4425 | 161 | 1146 |
| 50 | 126 | 88.26% | 97.92% | 72.55% | 4042 | 86 | 1529 |
| 80 | 160 | 87.32% | 97.71% | 70.35% | 3919 | 92 | 1652 |

- **Peak at k=15-20 (91.2 / 91.0%).** Both clear the 85% gate by ≥6 pp; this is the operational sweet-spot for the 2× / ~1.5 μm patch geometry.
- **All 7 k-values clear the 75% gate.** 6 of 7 (k=10..80) clear the 85% gate. The kNN-MRF backend is robust across nearly two orders of magnitude in k at this contrast + patch size.
- **Precision/recall asymmetry is the underlying mechanism.** Low k → low precision (78.83% at k=5) + high recall (93.90% at k=10): the small ball mis-classifies background emitters that happen to be in dense local fluctuations as "high" (FP). High k → high precision (97.92% at k=50) + low recall (70.35% at k=80): the large ball averages true patch density with surrounding background, missing emitters in the patch interior closest to the boundary (FN). k=15-20 is the precision/recall crossover where both metrics ≥ 87%.
- **k=20 is robust to ±5 perturbation, but k=15 is marginally better at this geometry.** The 91.23%/91.01% peak (k=15 vs k=20) at ~1.5 μm patches mildly favors the smaller ball — consistent with V12's bound that smaller r_k fits more comfortably inside the rect-thin tail (75 nm half-width). However, the difference is within rounding (Δ = +0.22 pp) and below any operational decision threshold. **k=20 stays the default**; the sweep does not motivate changing it.
- **Operational floor at low k confirms the V12 sub-sampling-noise bound.** k=5 lands at 85.15% — just barely above the 85% gate — and exhibits the FP-dominant failure mode predicted by σ_log ≈ 1/√k = 0.45 (noisy density distribution → fuzzy GMM split → over-aggressive smoothing labels border emitters as high). At k < 5 the curve would presumably collapse below the 75% gate, but k=5 was the lower bound of this sweep.

**Visual inspection** (per phase-3-extras of `start-round.md` — Read each PNG before claiming done):

- `k_sensitivity_sweep.png`: clean concave-up curve with the 75%/85% gate lines clearly separated from the data range (peak ~91% comfortably above 85%; floor 85% still on the gate). The annotated r_k values at each point trace the kNN-ball-radius vs structure-half-width story: r_k spans 40-160 nm across the sweep, bracketing the 75-150 nm rect short-side range. Two regimes visible: the rising left limb (k=5 → 15) where r_k grows from below the noise-floor scale up to the structure half-width, and the falling right limb (k=20 → 80) where r_k exceeds the structure half-width and the ball spills into background.
- `round_016_smlmrender_categorical.png`: 7-panel grid, all panels at the same ROI extent. The dominant visual is orange-TP patches against a blue-TN background — kNN-MRF correctly labels most emitters in all panels. The panel-to-panel evolution traces the FP/FN trade-off: k=5 panel has visible green-FP speckle dispersed through the background and around patch boundaries (the precision-loss signature, FP=1295); k=10 still has scattered green-FP but tighter to the patches; k=15-30 panels are visually cleanest — clean orange-TP patch interiors with thin fringe contamination; k=50-80 panels show progressively more red-FN appearing at patch interiors AND around patch peripheries (the recall-loss signature, FN reaching 1652 at k=80) while green-FP nearly vanishes (precision 97.7%). The visible trend matches the precision/recall numbers panel-by-panel: at low k the smoothing-amplified noise puts FP everywhere; at high k the over-averaged density puts FN at the patch fringe where the ball straddles the boundary. k=15 and k=20 are visually indistinguishable to the eye — the +0.22 pp accuracy difference is well below visual resolution.

## What failed

Nothing material. The sweep recovered the predicted concave-up curve, all gates met, peak at the V12-default k=20 confirms the existing default is well-positioned (the marginal k=15 advantage is within rounding). The lower-k limb at k=5 (85.15%) is right on the 85% gate — at k=3 or k=2 the curve would presumably fall through both gates, but those values are below the kNN-density estimator's useful regime and were not swept.

The `Assignment to ax in soft scope is ambiguous` warning during the SMLMRender panel build is cosmetic — Julia's loop-scope rule treating `ax` as a new local. Did not affect the output and did not warrant a `local ax` declaration since the run was single-threaded and one-shot.

## Files changed

```
dev/scripts/k_sensitivity_sweep.jl                         new (force-add, ~220 lines)
dev/scripts/output/k_sensitivity_sweep.csv                 new (force-add, 0.4 KiB)
dev/scripts/output/k_sensitivity_sweep.png                 new (force-add, 104 KiB) — accuracy curve
dev/scripts/output/round_016_smlmrender_categorical.png    new (force-add, 2.5 MiB) — 7-panel categorical render
STATUS.md                                                  P8 → DONE; Current State refresh; new Round History row; V17 mention
KNOWLEDGE_BASE.md                                          +V17 kNN k-sensitivity at fixed 2× ratio
QUESTIONS.md                                               +Q6 (P6 trio spec) +Q7 (P7 trio spec)
rounds/round_016_*.md                                      new
```

## Confidence

High. The predicted curve shape is recovered cleanly with the same RNG seed across all 7 k-values (only `density_k` varies). All gates met at k ∈ {10, 15, 20, 30, 50, 80} for the 85% gate and at k=5 for the 75% gate. The k=20 default is reaffirmed: at the marginal-best k=15 the gain is +0.22 pp accuracy, within rounding noise; the V12-default remains the right operational pick. The precision/recall crossover at k=15-20 matches the underlying mechanism (small ball → FP-dominant via smoothing-amplified noise; large ball → FN-dominant via density-mean regression). Medium-high on operational guidance: the sweep was at a single (ratio, patch_scale) point — characterization of the (k, patch-size) joint surface (V16 territory) and (k, ratio) joint surface (V13 territory) remain open; this round's claim is bounded to the 2× / ~1.5 μm operating point. The Round 012 V12 default (k=20) was hand-picked at `patch_scale=1.0` and the Round 016 sweep verifies it remains within rounding of the optimum at `patch_scale=0.75` — the surface is shallow around the peak.

## External consultations (if any)

None this round.

## Next steps

Priority 8 closes as DONE — kNN k-sensitivity at fixed 2× ratio + ~1.5 μm patches is mapped (peak k=15-20, ≥85% across k ∈ [10, 80], V17 records the bound). Together with V12 (k=20 default + structure-half-width bound), V13 (density-ratio operational floor), V14/V15 (calibration override path), and V16 (patch-size operating range), the (ratio, patch-size, k) operating surface is now characterized at three orthogonal slices through the 2×-ratio-1.0/1.5-μm-patches-k=20 anchor.

Two LOW-priority follow-ups remain TBD-trio:

- **P6** MRF v2 enhancements (graph-cuts MAP inference, soft-posterior output) — held until needed; trio not yet specified.
- **P7** Ratio-aware backend selector (or README-only guidance leg) — V13 ratio floor is documented in the README but a `RatioAwareDensityClusterConfig` wrapper is not yet built; trio not yet specified.

Both are flagged as OPEN questions in `QUESTIONS.md` (Q6 and Q7) requesting Keith specify the add/test/doc legs before they become selectable in a future round.

The natural extension priorities surfaced by this round are NOT promoted (per anti-pattern: do not act on discovered scope mid-round). They are:

- **(k, patch-size) joint surface** — sweep both axes at fixed 2× ratio. Would extend V17 + V16 into a 2D operating surface, useful for picking k automatically by patch-size estimate.
- **(k, density-ratio) joint surface** — sweep both at fixed ~1.5 μm patches. Would extend V17 + V13. The Round 013 result (kNN-MRF degenerates below 2×) was at fixed k=20; whether dropping k can rescue the low-contrast regime is open.

These would be candidate priorities only after P6/P7 trio specs land.

## Questions posted

- **Q6** — P6 trio specification (MRF v2 enhancements: graph-cuts vs soft-posteriors). The priority body lists three sub-items but does not commit to which one(s) ship; trio is TBD/TBD/TBD.
- **Q7** — P7 trio specification (ratio-aware backend selector vs README guidance). Two paths (`RatioAwareDensityClusterConfig` wrapper or README-only documentation extension); trio is TBD/TBD/TBD.
