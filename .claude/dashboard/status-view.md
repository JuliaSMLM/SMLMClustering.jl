# SMLMClustering

## Project Overview

Julia package providing a common interface over four single-molecule localization clustering backends (density-based, hierarchical density-based, Voronoi tessellation, and agglomerative hierarchical). Operates on SMLMData emitter tables and writes per-localization cluster labels back into the emitter struct so downstream rendering and analysis can treat cluster membership as a first-class attribute. Built to slot into the SMLMAnalysis pipeline after drift correction / grouping and before render.

## Development Overview

Fresh-scaffolded repo — only the `PkgTemplates` skeleton exists; no source yet. The interface shape has been agreed with the SMLMAnalysis coordinator: a single `cluster` entry point returning `(smld, ClusterInfo)`, shared config fields common to every backend, and SMLMData as the only upstream dependency. Next rounds build this out backend-by-backend starting with the abstract types and the density-based method.

## Task List

### Recent

| N | Focus | Status |
|---|-------|--------|
| 000 | Initial scaffold | done |
| — | — | — |
| — | — | — |

### Last round

The fresh-context round system was installed on top of the `PkgTemplates` skeleton, and the package's interface was pinned down through a design conversation with the SMLMAnalysis coordinator before any code got written. Durable knowledge about the label convention, the entry-point signature, the shared config layout, and the single-dependency rule now lives in the knowledge base so every future round starts from that agreement rather than rediscovering it.

**Technical detail:** Round-init scaffolded `.claude/commands/{start-round,dispatch-round}.md`, `.claude/scripts/dispatch-round.sh`, and `.claude/round/*.md` canonical refs + project-editable stubs; seeded `STATUS.md` / `KNOWLEDGE_BASE.md` / `QUESTIONS.md` at repo root. KB V1–V4 record: cluster labels in `emitter.id` (0 = noise, 1..K = cluster); `cluster(smld::BasicSMLD, cfg::AbstractClusterConfig) -> (smld, ClusterInfo)`; shared fields `min_points`, `use_3d`, `per_dataset`, `remove_unclustered` on every concrete config; Project.toml depends on SMLMData only. Phase 1 extras added for agent-inbox reads; Phase 4 extras send a one-line recap to `@analysis` after commit. Dashboard 2x2 tmux layout up; `/loop 3,33 * * * *` cron `08c747ed` dispatches `/dispatch-round` every 30 min.

<!-- ACTIVE -->

### Next

1. **[HIGH]** Scaffold package skeleton — add SMLMData dep, define `AbstractClusterConfig`, `ClusterInfo`, and the `cluster(smld, cfg) -> (smld, ClusterInfo)` entry point with no backends wired yet — TODO
2. **[HIGH]** Implement the density-based (DBSCAN) backend — shared fields plus `eps_nm`; dispatch writes labels into `emitter.id` with per-dataset handling — TODO
3. **[HIGH]** Implement the hierarchical density-based (HDBSCAN) backend — shared fields plus `min_cluster_size` — TODO
