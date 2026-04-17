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

**Decision:** Algorithm-specific distance fields on configs (`DBSCANConfig.eps_nm`, future `HierarchicalConfig.cut_nm`, etc.) are named with the `_nm` suffix and specified in nanometers. Emitter coordinates (`x`, `y`, `z`) on SMLMData emitter structs are in microns. Each backend converts internally via `radius_μm = eps_nm / 1000.0` before calling its underlying library.
**Evidence:** Round 002 `DBSCANConfig.eps_nm=100.0` correctly clusters at 100 nm on emitters with microscale coordinates. The `_nm` suffix makes the unit explicit at the call site so users don't accidentally pass microns. SMLMData's emitter docstrings state microns; Clustering.jl is unit-agnostic (takes a `Real` radius), so the conversion must live in our backend.
**Scope:** Every distance-valued config field names the unit in its suffix; backend code performs the conversion.

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
