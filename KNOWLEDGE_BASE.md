# KNOWLEDGE BASE — SMLMClustering

*Institutional memory across rounds. Validated Decisions accumulate; Dead Ends prevent repeating mistakes; Audits track KB health.*

*Read in Phase 1. Updated in Phase 4. Entries are never silently deleted — contradicted entries are marked RETIRED with back-references. See `.claude/round/review-protocol.md` for the KB audit procedure.*

---

## Validated Decisions

<!-- Durable claims with supporting evidence. Add as rounds validate design decisions. -->
<!-- Format:
### V1 — <short title>

**Decision:** <the claim>
**Evidence:** <what established this — round reference, test outcome, benchmark, etc.>
**Scope:** <when this applies; any limits>
-->

### V1 — Cluster labels live in `emitter.id`

**Decision:** Each backend writes cluster labels to `emitter.id` — `0` means noise/rejected, `1..K` are cluster ids. Optional `remove_unclustered` drops `id == 0` emitters from the output SMLD.
**Evidence:** Agreed with @analysis (SMLMAnalysis) at round-init; mirrors the free-integer-tag role `id` already has and avoids colliding with `track_id` (which has `0 = unlinked` semantics from FrameConnect).
**Scope:** All four backends (DBSCAN, HDBSCAN, Voronoi, hierarchical). Cluster step runs AFTER FrameConnect/BaGoL, so FrameConnect's use of `track_id` is not disturbed.

### V2 — Entry point is `cluster(smld, cfg) -> (smld, ClusterInfo)`

**Decision:** Single entry point `cluster(smld::BasicSMLD, cfg::AbstractClusterConfig) -> (smld, ClusterInfo)`. Dispatch on the concrete config type selects the backend.
**Evidence:** Tuple-return matches the ecosystem convention used by SMLMBoxer / BaGoL / Drift / FrameConnect / Render. @analysis wraps SMLMClustering via `const DBSCANConfig = SMLMClustering.DBSCANConfig` (upstream-owns-config rule), so the SMLMAnalysis `analyze()` dispatch re-exports without reimplementing.
**Scope:** Applies to all backends. `AbstractClusterConfig <: SMLMData.AbstractSMLMConfig`; `ClusterInfo <: AbstractSMLMInfo`.

### V3 — Shared config fields live on `AbstractClusterConfig`

**Decision:** Fields common to all backends — `min_points`, `use_3d`, `per_dataset`, `remove_unclustered` — are defined on every concrete `*Config` struct. Algorithm-specific fields (e.g. `eps_nm`, `cut_nm`, `density_factor`, `min_cluster_size`) live only on the struct that needs them.
**Evidence:** Agreed with @analysis at round-init. `per_dataset=true` by default — multi-dataset SMLDs get clustered within each dataset so `(dataset, id)` is unique.
**Scope:** All four concrete config structs.

### V4 — Upstream dep is SMLMData only

**Decision:** SMLMClustering depends only on SMLMData (plus whatever clustering libraries each backend needs, e.g. Clustering.jl). It does NOT depend on SMLMAnalysis — SMLMAnalysis is a downstream consumer that re-exports the configs in `src/steps/cluster.jl`.
**Evidence:** Matches SMLMBoxer / SMLMFrameConnection dep pattern. Keeps SMLMClustering usable as a pure algorithm package.
**Scope:** Project.toml `[deps]` must not include SMLMAnalysis.

### V5 — Config structs use `Base.@kwdef` with explicit defaults

**Decision:** Concrete `*Config` structs are declared with `Base.@kwdef struct Foo <: AbstractClusterConfig ... end`. Required-by-user fields (e.g. `DBSCANConfig.eps_nm`) have no default; shared fields use the defaults agreed in V3 (`min_points=5`, `use_3d=false`, `per_dataset=true`, `remove_unclustered=false`).
**Evidence:** Round 002 implemented `DBSCANConfig` this way. Keyword-only construction avoids positional-argument confusion as the shared-fields list grows across four backends, and the pattern composes cleanly with SMLMAnalysis's `const DBSCANConfig = SMLMClustering.DBSCANConfig` re-export (callers use `DBSCANConfig(eps_nm=50.0)` either way).
**Scope:** All four backend config structs. The `use_3d`/`per_dataset`/`remove_unclustered` defaults should be identical across backends so callers can swap backends without surprise.

### V6 — Distance units: configs take **nm**, emitter coords are **μm**

**Decision:** Algorithm-specific **distance-valued** config fields (`DBSCANConfig.eps_nm`, etc.) are named with the `_nm` suffix and specified in nanometers. Emitter coordinates (`x`, `y`, `z`) on SMLMData emitter structs are in microns. Each backend converts internally via `radius_μm = eps_nm / 1000.0` before calling its underlying library.
**Evidence:** Round 002 `DBSCANConfig.eps_nm=100.0` correctly clusters at 100 nm on emitters with microscale coordinates. The `_nm` suffix makes the unit explicit at the call site so users don't accidentally pass microns. SMLMData's emitter docstrings state microns; Clustering.jl is unit-agnostic (takes a `Real` radius), so the conversion must live in our backend.
**Scope:** Every unambiguously-distance-valued config field names the unit in its suffix; backend code performs the conversion. **Caveat (Round 011, Q5):** `HierarchicalConfig.cut_threshold` deliberately drops the `_nm` suffix because its unit depends on `linkage`: for distance-based linkages (`:single`, `:complete`, `:average`) it IS in nm and the backend divides by 1000; for `:ward` the dendrogram height is a variance-increase cost (roughly μm²) and is passed through unchanged. Fields whose unit is linkage- or algorithm-dependent must not carry the `_nm` suffix — pair them with a docstring that spells out the per-linkage interpretation, or provide an alternative parameter (here, `n_clusters`) that avoids the unit question.

### V7 — Voronoi backend is 2D-only; `use_3d=true` raises `ArgumentError`

**Decision:** `VoronoiConfig` does not implement 3D clustering. `cluster(smld, cfg::VoronoiConfig)` throws `ArgumentError` when `cfg.use_3d == true`, with a message directing users to `DBSCANConfig` or `HierarchicalConfig` for 3D data. The other three backends' `use_3d=true` paths remain supported.
**Evidence:** Round 004 adopted `DelaunayTriangulation.jl` for Voronoi tessellation; its 3D tessellation is not generally available, and every other pure-Julia Voronoi library (VoronoiCells.jl, VoronoiDelaunay.jl) is 2D only. Vendoring or building a 3D Voronoi from scratch is out of scope for a lightweight dep footprint (V4). Loud error prevents silent 2D fallback on 3D data.
**Scope:** Applies only to `VoronoiConfig`; the shared `use_3d::Bool = false` field stays on all four configs for uniformity, but the Voronoi dispatch rejects `true`.

### V8 — Voronoi density test: "dense" ⇔ cell area < mean / `density_factor`

**Decision:** `VoronoiConfig` defines a localization as dense when its Voronoi cell area is **strictly less than** `mean_cell_area / density_factor`, where the mean is computed over the clipped tessellation within the group (per `per_dataset`). Dense localizations connected via the Delaunay adjacency graph form raw clusters; clusters with fewer than `min_points` members are relabeled noise and the remainder are renumbered compactly `1..K` within the group.
**Evidence:** Round 004: three tight σ=10 nm blobs amid 60 scattered noise points are recovered as three clusters at `density_factor=2.0`; lowering `density_factor` below 1 inflates the threshold and absorbs most noise; raising `density_factor` shrinks the threshold and demands tighter packing. Matches the SR-Tesseler formulation (Levet et al., Nat Methods 2015). Polygons are clipped to the convex hull (`voronoi(tri; clip=true)`) so every generator has a finite area, with a documented caveat that hull cells are systematically smaller than their infinite-plane area.
**Scope:** `VoronoiConfig` only. `density_factor > 0` is validated; `density_factor → 0` makes the threshold diverge (everything is "dense") and one giant component can form, while very large `density_factor` leaves few or no dense points.

### V9 — `cluster()` is non-mutating — input SMLD is never modified

**Decision:** Every backend's `cluster(smld, cfg)` method begins by deep-copying the input emitter vector and wrapping it in a fresh `BasicSMLD` (sharing the input's camera, frame/dataset counts, and metadata references). All label writes (`emitter.id = ...`) and downstream filtering (`remove_unclustered`) operate on the copy. The returned SMLD is the one backends wrote into; the caller's original SMLD is unmodified.
**Evidence:** Round 011 processed `QUESTIONS.md` Q3 (answered by Keith 2026-04-17: "we don't modify our input; that is the interface convention used in SMLMAnalysis"). Downstream SMLMAnalysis steps (detectfit, filter, frameconnect, drift, render, densityfilter, intensityfilter) universally return new BasicSMLDs rather than mutating input; non-mutating `cluster()` aligns with that ecosystem convention. Implementation choice `deepcopy(smld.emitters)` (per @analysis's option (b)) is backend-agnostic — any `<: AbstractEmitter` subtype works without per-type `_with_dataset`-style constructors. Non-mutation is asserted in every backend's tests: after `cluster(smld, cfg)`, the input `smld` still has all emitter `id == 0`.
**Scope:** Applies to all current and future `cluster(smld, cfg::AbstractClusterConfig)` methods. A hypothetical `cluster!` mutating variant would be a separate method with the `!` suffix; the un-suffixed entry point must never mutate. Note this contrasts with the `cluster_statistics(smld, cfg)` entry point (V10), which is a passthrough (SMLD reference is the same as input) because it writes nothing — that asymmetry is intentional and is documented on both docstrings.

### V10 — `cluster_statistics` sibling interface, pass-through SMLD, summary-scalar + extras convention

**Decision:** A second entry point `cluster_statistics(smld, cfg::AbstractStatisticsConfig) -> (smld, ClusterStatisticsInfo)` parallels `cluster()` for read-only spatial-statistic backends (clustering tendency, density diagnostics, ...). Three structural choices that govern all current and future backends in this hierarchy:

1. **Pass-through SMLD reference (NOT a deep-copy):** `cluster_statistics` writes nothing onto the SMLD and returns the *same reference* as the input. The two-tuple return shape is preserved for ecosystem symmetry, but no allocation or mutation occurs. This is intentionally asymmetric with `cluster()`'s V9 deep-copy guarantee — copying for a read-only operation is wasted work.
2. **Abstract supertype + dispatch on concrete config:** `AbstractStatisticsConfig <: SMLMData.AbstractSMLMConfig` is the abstract supertype; concrete backends (e.g. `HopkinsConfig`) supply a `cluster_statistics(smld, cfg::SomeConfig)` method. The fallback errors with a message naming the available concrete backends.
3. **Summary scalar + extras dictionary:** `ClusterStatisticsInfo` carries a single `statistic::Float64` (the headline number — `info.statistic` ergonomic for one-number consumers) plus `extras::Dict{Symbol,Any}` for vector-valued or supplementary outputs. Backends that produce a natural vector output put the vector in `extras` under a descriptive key (e.g. `:hopkins_per_dataset`, `:density_per_emitter`) and a meaningful aggregate (mean, median) in `statistic`. This avoids each backend reinventing its own info-struct shape.

Shared fields on every concrete `*StatisticsConfig` struct: `use_3d::Bool=false`, `per_dataset::Bool=true`. Algorithm-specific fields live on the concrete subtype.

**Evidence:** Task A direct execution 2026-04-27 landed P11+P12+P13+P14 as a single bundle: `AbstractStatisticsConfig`, `ClusterStatisticsInfo`, `cluster_statistics` fallback, and `HopkinsConfig` backend. Test suite verifies the pass-through guarantee (`smld_out === smld` after a `cluster_statistics` call), the summary-scalar+extras convention (Hopkins per-dataset vector lives in `extras[:hopkins_per_dataset]`, mean across datasets in `statistic`), and the asymmetry with `cluster()` (which still deep-copies). 221/221 tests pass.

**Scope:** Applies to all current and future `cluster_statistics(smld, cfg::AbstractStatisticsConfig)` methods. The pass-through guarantee is binding — backends must not mutate the SMLD or its emitters. Future statistic backends (Voronoi density, Ripley K, NND distribution, ...) all subtype `AbstractStatisticsConfig` and follow the summary-scalar+extras convention. NaN is the canonical "couldn't compute on this group" return (empty SMLD, n_samples > n_points, degenerate bbox) — backends never throw on data-shape edge cases that are valid but produce no statistic; argument validation (n_samples ≥ 1, etc.) still throws `ArgumentError` at the boundary.

### V11 — MRF density-regime clustering pipeline (`MRFDensityClusterConfig`)

**Decision:** Adaptive-density clustering pipeline for data with multiple density regimes is structured as four serial steps:

1. **Per-emitter Voronoi density** → log ρᵢ = log(1/Aᵢ). Reuses the helper `_voronoi_areas` in `src/utils.jl` (extracted from `voronoi_density.jl` so both `VoronoiDensityConfig` and `MRFDensityClusterConfig` share the tessellation logic).
2. **Regime assignment.** Either an explicit `regime_thresholds::Vector{Float64}` of length `n_regimes - 1` (binning, GMM bypassed; per-emitter unary `U[i, k]` is `0` for the matching bin and `1e6` for others), OR an `n_regimes`-component 1D Gaussian mixture EM on log ρ producing unaries `U[i, k] = -log(w_k * N(log_rho[i] | μ_k, σ_k²))`. GMM components are sorted ascending by mean — convention is **regime 1 = lowest density (treated as background/noise) and `n_regimes` = highest density**.
3. **Multi-class Potts MRF refinement via ICM** over the Delaunay neighbor graph (default, free since step 1 already triangulated) or a symmetrized k-NN graph. Pairwise term `V(xᵢ, xⱼ) = 0 if xᵢ = xⱼ else 1`. Auto-tuned smoothness `λ = max(1e-6, MAD(U_max_i - U_min_i))` per group when `smoothness_lambda === nothing`. ICM iterates until no point changes label or until `icm_iters` is reached. Only `:icm` is supported in v1; `:graph_cuts` is the future-extension slot but raises `ArgumentError` today.
4. **Connected components** via BFS on the same neighbor graph restricted to foreground nodes (regime ≥ 2). Components below `min_points` are demoted to noise (id = 0).

Per-cluster outputs land in `emitter.id` (V1 convention). MRF-specific outputs go on `smld_out.metadata` (mirrors HDBSCAN's metadata-stamping pattern from `hdbscan.jl`):
- `metadata["mrf_regime_per_emitter"]::Vector{Int}` — per-emitter regime in 0..`n_regimes`, in original emitter order (0 = ungroupable / group too small / GMM failed).
- `metadata["mrf_lambda_used"]::Vector{Float64}` — per-group λ actually used, in `_group_by_dataset` order.
- `metadata["mrf_regime_means"]::Vector{Vector{Float64}}` — per-group GMM means (sorted ascending), filled with `NaN`s when `regime_thresholds` was supplied (signals "manual binning, no fit").

**Evidence:** Direct execution 2026-04-29 landed `MRFDensityClusterConfig` per Keith's design picks (hard regime labels, multi-regime from start, GMM, per-dataset with override). 981/981 tests pass with `SMLM_TEST_FULL=1` (39 new MRF tests covering 2-regime / 3-regime / missing-middle / spurious-small / threshold-override / per_dataset / :knn / determinism / edge cases). Fast tier 161/161 in 33s. Validates the four-step pipeline against synthetic data exhibiting both target failure modes (false-positive small knots in a low-density sea + missing middles in genuine dense regions).

**Scope:** Applies to `MRFDensityClusterConfig` and any future density-regime variants (e.g. graph-cuts inference replacing ICM). 2D only (V7 — DelaunayTriangulation.jl). Lowest-regime-as-noise convention is binding — backends that subclass this pipeline must not invert it without explicit user-facing rename. The auto-λ heuristic (MAD of unary range) is data-scale-free and works without dataset-specific tuning; users with strong priors override via `smoothness_lambda`.

### V13 — kNN-MRF operational floor on density ratio: clean at ≥ 2×, degenerates below; voronoi-GMM is the better backend at low contrast

**Decision:** On a controlled 5×5 μm A431-mimic synthetic with fixed patch geometry, the `MRFDensityClusterConfig` kNN-density pipeline (`density_estimator=:knn, density_k=20, n_regimes=2`) has a sharply asymmetric operational regime relative to the high-density / low-density density ratio:

- ratio ≥ 2.0× — kNN-MRF dominates all alternatives (DBSCAN, HDBSCAN, voronoi-GMM); accuracy 89-96% across {2, 3, 5}×.
- 1.65× ≤ ratio < 2.0× — kNN-MRF clears the 75% accuracy gate but is roughly tied with voronoi-GMM.
- ratio < 1.65× — kNN-MRF degenerates: at ratio 1.2× it collapses to 35.4% accuracy by labeling nearly the entire ROI "high" (FP catastrophe, not FN). Voronoi-GMM (no MRF smoothness) holds 64.5% at ratio 1.2× and is the recommended backend for low-contrast datasets.

The 85% gate falls at ratio ≈ 1.85×; the 75% gate at ratio ≈ 1.65×, both by linear interpolation between the swept ratio points.

**Evidence:** Round 013 sweep at ratios {1.2, 1.5, 2.0, 3.0, 5.0}× — script `dev/scripts/density_ratio_sweep.jl`, results CSV at `dev/scripts/output/density_ratio_sweep.csv`, aggregate curve at `dev/scripts/output/density_ratio_sweep.png`, per-backend categorical TP/TN/FP/FN circle-plot grids. Same patch geometry across ratios (controlled `SIM_SEED`); only `RHO_HIGH_BONUS` varies. Sanity check: 2.0× column reproduces Round 012's 89.11% headline exactly (same patches, same RNG, same backend) — refactor of `simulate_a431_mimic.jl` to expose `simulate_dataset(; rho_low, rho_high_bonus, seed, verbose)` does not change semantics. DBSCAN trails throughout (always 100% recall, terrible precision: at min_points=5/eps=100nm it labels nearly everything "high"). HDBSCAN flat ~50-65% (its parameter regime doesn't engage with this density-segmentation task). Visual `density_ratio_sweep_mrf_knn_panels.png` row confirms the 1.2× failure mode is FP saturation across the whole ROI, not interior FN.

The 1.2× failure mechanism: when high-density and low-density are close, the kNN log-density distribution has two modes that nearly overlap, the GMM EM fit picks a degenerate split (one mode swallows most of the data), and the multi-class Potts smoothness term then propagates that majority label across neighbor edges via ICM until the entire foreground graph is labeled as one class. Voronoi-GMM avoids this because it is per-emitter independent — without the smoothness prior, individual emitters get classified on their own log-ρ value and the failure is graceful (low recall) rather than catastrophic (FP saturation).

**Scope:** Applies to `MRFDensityClusterConfig` with the 2-regime, kNN-density configuration on patch-and-background data resembling the A431-mimic geometry (1-3 μm patches, AR 1-20, 5×5 μm ROI, RHO_LOW = 500/μm²). For real dSTORM with similar density character, deploy MRF when the rough density ratio (estimable from the data via a histogram of voronoi-density log-ρ) exceeds 2×; below that, prefer voronoi-GMM or another non-smoothing density classifier. Higher `n_regimes` and lower `density_k` may shift the operational floor — those are not characterized here. The Potts-overpowers-GMM-signal mechanism is intrinsic to MRF + EM-fit unaries; alternative inference (graph cuts, soft posteriors) might rescue the low-contrast regime but is out of scope for V13.

### V14 — `regime_thresholds` override is the right tool at low contrast (where GMM auto-mode degenerates), but forfeits the Potts-smoothness benefit at high contrast

**Decision:** The undocumented `MRFDensityClusterConfig.regime_thresholds::Vector{Float64}` override hard-pins the unary cost matrix to `0` for the matching bin and `1e6` for others (`_unary_from_thresholds` in `src/backends/mrf_density.jl`). On the synthetic A431-mimic, this produces an asymmetric operational regime relative to the high-density / low-density density ratio:

- ratio ≥ 2× — **GMM auto-mode wins by ~5 pp** (88.2% override vs 93.7% auto on the held-out 2× query ROI). The hard `1e6` penalty dwarfs the Potts pairwise term `λ × #neighbors_with_label_not_k`, so a borderline interior emitter whose individual log ρ falls below the threshold by chance can no longer be rescued by its high-density neighbors. The Potts-smoothness benefit that the kNN-MRF default depends on (V12 patch-interior FN-band collapse) is forfeit; the +784 increase in FN at almost equal FP rate is the smoothness-loss signature.
- ratio < 2× — **override wins by ~4 pp** (75.9% override vs 72.1% auto on the 1.5× query ROI calibrated on a 2× ROI). At low contrast GMM auto-mode degenerates (the two log ρ modes coalesce, EM picks a near-tied split) and the soft Potts smoothness amplifies the degeneracy into uniform misclassification (V13 mechanism). Hard binning prevents the amplification: even though the calibrated threshold is "wrong" for 1.5× data (naturally too high), it produces a graceful drop in recall (33.3% vs 52.6%) rather than the catastrophic FP saturation auto-mode falls into.

**Operational guidance.** Use override mode when (a) calibration ROIs exist with stable density structure across query ROIs, AND (b) data is low-contrast (ratio < 2×) where GMM auto-mode is at risk of degenerate fits. Use GMM auto-mode at high contrast (ratio ≥ 2×) where the Potts smoothness term materially helps. The 2× ratio crossover is empirical and may shift with `density_k`, `n_regimes`, or different patch geometry; the underlying mechanism (hard penalty dwarfs pairwise term) is stable. The clean fix is a soft-unary override mode that lets users supply `(means, vars, weights)` for each regime — producing the same soft Gaussian unaries that GMM auto-mode produces, preserving the Potts smoothness benefit AND skipping per-ROI EM (Pareto improvement; logged as STATUS Priority 9).

**Helper API.** `calibrate_regime_thresholds(smld; n_regimes=2, density_estimator=:knn, density_k=20)::Vector{Float64}` exported from `src/backends/mrf_density.jl`. Computes per-emitter log density via Voronoi or kNN, fits 1D `n_regimes`-component GMM via `_gmm_em_1d`, returns analytic Bayes decision boundaries between consecutive (mean-sorted) components — the quadratic root in `[μ_k, μ_{k+1}]` of `log(w_k) + log N(x | μ_k, σ_k²) = log(w_{k+1}) + log N(x | μ_{k+1}, σ_{k+1}²)`, with midpoint fallback when the discriminant is negative. Output shape is exactly what `MRFDensityClusterConfig.regime_thresholds` accepts.

**Evidence:** Round 014 calibration demo at `dev/scripts/regime_thresholds_calibration.jl` — calibration ROI seed 20260429 (2× ratio), query ROIs seed 20260430 (2× and 1.5× ratios), single threshold fit at log ρ = 6.685. Outputs `dev/scripts/output/regime_calibration_demo.{csv,png}` and `dev/scripts/output/regime_calibration_panels.png` (4-panel categorical TP/FP/FN circle plot, one per case). Visual diagnostic confirms the mechanism panel-by-panel: 2× override panel has visibly denser cyan-FN at patch interiors (smoothness loss); 1.5× override panel has the orange-FP saturation collapsed and replaced by cyan-FN at patch interiors (Potts amplification avoided). Thorough-tier tests cover argument validation, 2-regime + 3-regime fits, kNN path, and end-to-end use of returned thresholds in `MRFDensityClusterConfig`. Test suite green at fast 172/172 / thorough 1002/1002.

**Scope:** Applies to `MRFDensityClusterConfig` 2-regime kNN-density configuration on patch-and-background data resembling the A431-mimic geometry. The `±2 pp at 2×` pre-empirical gate that Priority 3 was specified against was not met; the asymmetric outcome IS the finding. The crossover ratio between "auto wins" and "override wins" sits in [1.5×, 2.0×] but was not finely swept here — Priority 9 (soft-unary override) and Priority 4 (patch-size sweep) are the natural extensions. Hard-threshold override is NOT recommended at high contrast unless calibration cost dominates (override skips per-ROI EM); for those cases the soft-unary follow-up (P9) is the right path.

### V12 — kNN density estimator beats Voronoi at thin-elongated-patch interiors

**Decision:** `MRFDensityClusterConfig` exposes `density_estimator::Symbol = :voronoi` (default, backward compatible) with `:knn` as an alternative paired with `density_k::Int = 20`. The kNN estimator computes ρᵢ = k / (π · r_kᵢ²) where r_kᵢ is the kth nearest-neighbor distance for emitter i. For SMLM datasets containing thin elongated structures (widths comparable to the local nearest-neighbor distance — e.g. AR ≥ 5 fibers narrower than ~150 nm at high density) the kNN estimator dramatically reduces the patch-interior false-negative band that Voronoi produces.

**Evidence:** Round 012 on the synthetic A431-mimic (`dev/scripts/output/synthetic_smld.jld2`, 5×5 μm ROI, 12 patches AR 1.1-19.7, 13,579 emitters). Voronoi MRF: 79.85% headline accuracy, 22.49% interior FN-rate (FN concentrated >100 nm inside patches). kNN MRF (k=20): 89.11% accuracy (+9.26 pp), 6.45% interior FN-rate (-16.04 pp). Per-patch breakdown shows the win is consistent across all patches with measurable interior emitters: rect AR 8.9 went from 35% → 2.5% interior FN; ellipses dropped 60-93%. `dev/scripts/output/mrf_interior_diagnosis.png` row-3 shows the FN signed-distance distribution shifted from peaking at -0.2 μm (deep interior) to peaking at 0 μm (boundary), with the deep-interior tail collapsed. Density distributions in row 1 are visibly more separable under kNN — σ_log shrinks roughly 1/√k as expected.

**Scope:** Applies to MRFDensityClusterConfig only. Operational regime bound: kNN density requires the kNN ball radius (≈ √(k / π · ρ)) to be smaller than the structure half-width, otherwise the ball spills into background and the GMM regime split can flip foreground/background. For thin patches where k=20 is too coarse, drop to k=8-10. The Voronoi default is preserved because Voronoi remains the right estimator for blob-shaped clusters where boundary spillage is a non-issue and backward compatibility matters. When `graph_kind=:delaunay`, the Voronoi tessellation is still computed (needed for the neighbor graph) regardless of the density estimator selection.

---

## Dead Ends

<!-- Approaches confirmed not to work, so future rounds don't repeat them. -->
<!-- Format:
### D1 — <short title>

**What was tried:** <approaches>
**Why it failed:** <concrete reason>
**Round reference:** round_NNN_slug.md
**Workaround (if any):** <alternative>
-->

### D1 — No lightweight registered Julia HDBSCAN library (as of 2026-04-17)

**What was tried:** `Clustering.jl` 0.15.8 (already a dependency) — no `hdbscan` symbol. Julia General registry search — no registered HDBSCAN package. `baggepinnen/HDBSCAN.jl` (GitHub-only) — PyCall wrapper requiring Python `hdbscan`; Python not available in this environment. `HorseML.jl` (registered, 0.4.1) — has a pure Julia HDBSCAN implementation, but the package pulls in CUDA, NNlib, NNlibCUDA, and Zygote; adding it as a dep violates the lightweight algorithm-package constraint (V4).
**Why it failed:** No registered Julia HDBSCAN package with acceptable dependency weight exists. `Clustering.jl` issue #139 tracks adding HDBSCAN but it has not merged as of v0.15.8.
**Round reference:** round_003_hierarchical-backend.md
**Workaround (if any):** (a) Wait for `Clustering.jl` to merge HDBSCAN (monitor issue #139). (b) Implement HDBSCAN from scratch using KNN + MST + hierarchy extraction (~300 lines, non-trivial). (c) Extract just the HDBSCAN source from `HorseML.jl` into a vendored file (check license). Option (a) is lowest risk; check again when `Clustering.jl` bumps past 0.15.8.

---

## Audits

<!-- Written by rounds that run the KB audit (every 2 × review_cycle_rounds). See review-protocol.md. -->

(none yet)
