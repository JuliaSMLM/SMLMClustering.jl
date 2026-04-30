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
