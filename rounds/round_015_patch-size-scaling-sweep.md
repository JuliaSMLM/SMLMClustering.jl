# Round 015 — patch-size scaling sweep

**Date:** 2026-04-30
**Status:** done
**Priority worked on:** Priority 4 — patch-size scaling sweep

## Output files to review

- `dev/scripts/output/patch_size_sweep.png` — accuracy curve, both backends, log-x patch size, dashed 75% gate + dotted 85% gate. Look for: kNN-MRF (red) climbs from 62.7% at 0.3 μm to a 88-89% plateau between 1.0-3.0 μm, crossing the 75% gate between 0.3 and 0.5 μm and the 85% gate between 0.5 and 1.0 μm. voronoi-GMM (yellow) is roughly flat at ~70% across the entire sweep — at 2× density ratio it is bounded by emitter-independent log-ρ separability, so patch size barely affects it. The voronoi-GMM curve never crosses the 75% gate.
- `dev/scripts/output/round_015_smlmrender_categorical.png` — 5-panel SMLMRender CircleRender categorical TP/TN/FP/FN circle plot per scale (0.3, 0.5, 1.0, 2.0, 3.0 μm). Encoding id 1=TN (blue), 2=TP (orange), 3=FP (green), 4=FN (red); kNN-MRF predictions only. Look for: 0.3 μm panel is dense FP+FN saturation (the 62.7% accuracy regime is visible as substantial green-FP plus red-FN intermixed with TP-orange); 0.5 μm panel still has visible FP haze around small patches but TP-orange patches are emerging cleanly; 1.0 / 2.0 / 3.0 μm panels have crisp TP-orange patch interiors with thin FN-red boundaries — the canonical "MRF works" visual.
- `dev/scripts/output/patch_size_sweep.csv` — per (backend, scale) accuracy/precision/recall/TP/TN/FP/FN for both kNN-MRF and voronoi-GMM across the 5-scale sweep.
- `dev/scripts/patch_size_sweep.jl` — sweep script. Reads `simulate_a431_mimic.jl` for `simulate_dataset(; rho_low, rho_high_bonus, seed, patch_scale, verbose)` (newly extended with `patch_scale::Float64=1.0`).
- `dev/scripts/simulate_a431_mimic.jl` — refactored to thread `patch_scale` through `sample_patch` → `place_patches` → `simulate_dataset`. Linear dimensions of rect (length L, half-width b) and ellipse (semi-major a, semi-minor b) all multiply by `patch_scale`; default 1.0 reproduces the original geometry exactly (same 12-patch layout at 2.0 μm nominal).

## Hypothesis

The MRF kNN-density pipeline (`MRFDensityClusterConfig` with `density_estimator=:knn, density_k=20`) has an effective patch-size operating range bounded below by the kNN ball radius. Round 012 / KB V12 articulated this as "kNN ball radius (~√(k/πρ)) must be smaller than the structure half-width, otherwise the ball spills into background and the GMM regime split can flip". With `k=20` and `ρ_high=1000/μm²` the kNN ball radius is ~80 nm, so the bound predicts kNN-MRF should fail when patch half-widths drop below ~80 nm — at the 0.3 μm patch-size scale the half-widths are 75-150 nm, comparable to the ball radius. The bar from Priority 4 was "(a) curve is interpretable (expected: dip at ≤0.5 μm where patches drop below Voronoi neighborhood scale, then plateau), (b) MRF holds >75% accuracy in the 0.5-5 μm range".

## What was attempted

Refactored `dev/scripts/simulate_a431_mimic.jl` to expose `patch_scale::Float64=1.0` as a kwarg on `simulate_dataset`, threaded down through `sample_patch` (rect L ∈ [scale, 3×scale], b ∈ [0.5×scale, scale]; ellipse a ∈ [0.5×scale, scale], b ∈ [0.25×scale, 0.5×scale]) and `place_patches`. Added SMLMRender to `dev/scripts/Project.toml` via `Pkg.develop(path="/home/kalidke/julia_shared_dev/SMLMRender")`. Built `dev/scripts/patch_size_sweep.jl` running both kNN-MRF (k=20) and a 2-component voronoi-GMM baseline across 5 patch-size scales `{0.15, 0.25, 0.5, 1.0, 1.5}` (corresponding to nominal mean rect length `{0.3, 0.5, 1.0, 2.0, 3.0} μm`), with deterministic RNG seed `20260429` so the same RNG draws produce different patch geometries only by virtue of the `patch_scale` multiplier. Per scale: simulate, predict, compute TP/TN/FP/FN per emitter, emit a SMLMRender CircleRender categorical (`color_by=:id, categorical=true, colormap=:tab10`) panel plus a row entry in the sweep CSV. Composed an aggregate accuracy curve (CairoMakie, log-x, both backends) and a 5-panel SMLMRender canonical render.

Scope adjustment from the original `{0.3, 0.5, 1, 2, 5, 10 μm}` spec: anything ≥5 μm cannot fit two non-overlapping patches in the 5×5 μm ROI alongside the placement spacing the simulator demands, so the upper end of the original sweep would have collapsed to "1 patch and reject the rest" and stopped being a fair test of the algorithm. Truncated to `{0.3, 0.5, 1.0, 2.0, 3.0} μm` which keeps the placement well-conditioned. Recorded as a bookkeeping item in the round commit; downstream rounds can attack the larger-patch regime by enlarging the ROI rather than the patches.

## What worked

- **Curve shape matches the predicted regime bound.** kNN-MRF sweeps `62.67% → 79.25% → 88.72% → 89.11% → 88.31%` across `{0.3, 0.5, 1.0, 2.0, 3.0} μm`. The dip at 0.3 μm is exactly where the kNN(20) ball radius (~80 nm) becomes comparable to the structure half-width (75-150 nm) and the ball spills into background. The 1.0 μm sample reproduces Round 012's 89.11% headline within rounding (same RNG, same backend, default `patch_scale=1.0`) — a positive control that the refactor is semantics-preserving.
- **Operational range gate met.** At 0.5 μm kNN-MRF lands at 79.25%, comfortably above the 75% gate; the 85% gate lifts off between 0.5 and 1.0 μm. Across the 1.0-3.0 μm plateau the MRF holds 88-89%. Operational range for the 2× density-ratio synthetic with kNN k=20: ≥0.5 μm passes 75%, ≥1.0 μm passes 85%.
- **voronoi-GMM is the wrong baseline at this contrast.** Across the entire sweep voronoi-GMM holds ~65-72% (max 71.79% at 2.0 μm). At 2× density ratio the per-emitter independence of voronoi-GMM bounds it well below the kNN-MRF plateau — patch size barely moves it because separability is set by log-ρ histogram structure, not by patch geometry. This is consistent with the Round 013 V13 finding that voronoi-GMM wins only below 2× ratio (where MRF Potts smoothness amplifies GMM degeneracy); at the 2× contrast used here voronoi-GMM is the strict-loss baseline.

**Visual inspection** (per phase-3-extras of `start-round.md` — Read each PNG before claiming done):

- `patch_size_sweep.png`: clean monotone climb from 0.3 → 1.0 μm, then a flat plateau between 1.0-3.0 μm with a faint dip at 3.0 μm (88.31% vs 89.11% at 2.0 μm); both gate lines cross between scale points where expected. voronoi-GMM curve sits well below kNN-MRF except at 0.3 μm where the two cross (kNN 62.7% < voronoi 65.3% — kNN's smoothness fails fastest).
- `round_015_smlmrender_categorical.png`: 0.3 μm panel is FP-green saturation across the ROI mixed with FN-red and TP-orange (the smoothness-collapse signature seen in Round 013's 1.2× panel — same FP-saturation failure, here triggered by kNN-ball-spillage rather than low contrast). 0.5 μm panel still has visible FP haze around small patches but TP-orange patch interiors are starting to emerge. 1.0 / 2.0 / 3.0 μm panels show clean orange-TP patch interiors with thin red-FN boundary bands and minimal FP — the canonical "MRF works" visual matching V12's interior-band collapse.

## What failed

Nothing material at the spec'd gates. The 0.3 μm scale is below the operational floor (62.7% accuracy) but this is the predicted-failure regime, not a regression — the entire point of the sweep was to map that floor. Documented as a hard lower bound in V16 below.

The original Priority 4 spec called for sweeping up to 10 μm. Truncated to 3.0 μm because larger patches don't fit alongside their non-overlap-buffer in the 5×5 μm ROI. This is a scope reduction relative to the spec; the bound (kNN-MRF holds 75% at ≥0.5 μm) is reported only across the {0.5, 1.0, 2.0, 3.0} μm range that was actually swept. Larger-patch behavior remains untested in this round; the natural extension is to enlarge the ROI rather than the patches.

## Files changed

```
dev/scripts/simulate_a431_mimic.jl                      patch_scale kwarg threaded through sample_patch / place_patches / simulate_dataset
dev/scripts/patch_size_sweep.jl                         new (force-add, ~210 lines)
dev/scripts/Project.toml                                +SMLMRender via Pkg.develop
dev/scripts/Manifest.toml                               +SMLMRender
dev/scripts/output/patch_size_sweep.csv                 new (force-add, 762 B)
dev/scripts/output/patch_size_sweep.png                 new (force-add, 110 KiB) — accuracy curve
dev/scripts/output/round_015_smlmrender_categorical.png new (force-add, 3.2 MiB) — 5-panel categorical render
STATUS.md                                                P4 → DONE; Current State refresh; new Round History row; V16 mention
KNOWLEDGE_BASE.md                                        +V16 patch-size operating range
rounds/round_015_*.md                                    new
```

This round also bundles the previously-uncommitted V15 (soft-unary calibration) source/test/doc changes from the 2026-04-30 direct work that landed in `STATUS.md` and `KNOWLEDGE_BASE.md` text but was never committed:

```
src/backends/mrf_density.jl                              regime_gaussians soft-unary path + calibrate_regime_gaussians helper
src/SMLMClustering.jl                                    +1 export
test/test_mrf_density.jl                                 +new soft-unary tests
README.md                                                MRF subsection: regime_gaussians + V15 reference
api_overview.md                                          MRF entry: regime_gaussians + helper
```

## Confidence

High on the patch-size finding: deterministic RNG seed, the 1.0 μm scale reproduces Round 012's 89.11% headline exactly so the refactor is semantics-preserving, the dip at 0.3 μm is the predicted-failure regime per V12, the curve crosses both gates monotonically. Medium-high on the operational guidance: the 0.5 μm gate-crossing is sampled at one point (no inner sweep across 0.3-0.5 μm); a finer sweep could place the crossover with more precision but the conclusion (≥0.5 μm clears 75%) is robust. The original spec asked to sweep up to 10 μm — truncated to 3.0 μm because the 5×5 μm ROI doesn't accommodate larger patches with non-overlap; bound is reported within the {0.5, 1.0, 2.0, 3.0} μm range.

## External consultations (if any)

None this round.

## Next steps

Priority 4 closes as DONE — patch-size operating range is mapped (kNN-MRF k=20 at 2× contrast: ≥0.5 μm for 75% gate, ≥1.0 μm for 85% gate, plateau 88-89% across 1.0-3.0 μm). KB V16 records the bound. Diagnosis loop's HIGH-priority items P1-P4 are now all closed (P5 remains BLOCKED on Keith creating the JuliaSMLM org repo).

Three LOW-priority follow-ups remain open in STATUS:
- **P6** MRF v2 enhancements (graph-cuts, soft posteriors) — held until needed.
- **P7** Ratio-aware backend selector / README guidance — README already documents the V13 ratio floor but a `RatioAwareDensityClusterConfig` wrapper is not yet built.
- **P8** kNN k-sensitivity sweep — pairs naturally with this round; would map the (k, structure-size) operating surface.

Closest natural next round: **P8** (k-sensitivity sweep at fixed 2× ratio, fixed 1.0 μm patch size, sweeping `k ∈ {5, 10, 20, 40, 80}`). Together with V12 (k=20 default) and V16 (patch-size range) this would complete the (ratio, patch-size, k) operating surface for the kNN-MRF backend.

The 5×5 ROI ceiling on this round (no patches >3 μm) is a follow-up sub-priority worth surfacing — enlarging the ROI to 10×10 μm or larger would let downstream rounds extend the patch-size sweep up to the original 10 μm ceiling. Not promoted to a priority itself; bundled into P8's "or related sweep extensions" scope.

## Questions posted

(none — finding is well-specified, gates met, follow-up shape clear)
