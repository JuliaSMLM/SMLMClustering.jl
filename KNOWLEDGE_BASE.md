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

(none yet)

---

## Audits

<!-- Written by rounds that run the KB audit (every 2 × review_cycle_rounds). See review-protocol.md. -->

(none yet)
