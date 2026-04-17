# QUESTIONS — SMLMClustering

*Async human-in-the-loop. Rounds post OPEN questions when they need human judgment (disputed review findings, dependency additions, design decisions). Human answers between rounds by editing items into ANSWERED. Next round processes ANSWERED and marks PROCESSED.*

---

## OPEN

<!-- Items awaiting human input. Each item MUST be written for a human who
     has not read any recent round files. Four-block shape, in this order:

### Q1 — <title that ends in a question mark?>

**From round:** NNN

**Short question:** <one plain-English sentence, answerable in isolation — no jargon, no variable names, no file paths>

**What each answer means:**
- **Yes / Option A:** <what changes in the project if the human picks this>
- **No / Option B:** <what changes if they pick this>
- (more options as needed)

**Technical detail:** <free-form context for future-round workers: variable names, file paths, flag values, reference commits, code snippets. This block is NOT read by the human answering the question; it exists so the round that PROCESSES the answer has enough context to act.>

HARD RULE: A human who has not read the round file should be able to answer the Short Question using only the What-each-answer-means block. If they need to read Technical detail to answer, the Short Question is wrong — rewrite it.
-->

### Q3 — Should the clustering entry point be renamed to signal that it modifies its input?

**From round:** 005

**Short question:** In Julia, functions that modify their arguments conventionally end in an exclamation mark (e.g. `sort!`). The current clustering entry point is named without one, even though it modifies the input data in place. Should it be renamed so users see at a glance that calling it changes their input?

**What each answer means:**
- **Rename it (with the `!`):** Users of the package see a clear signal in the name that the function modifies its input. Any code already using the old name has to be updated — today that is one downstream consumer with a one-line re-export, so the cost is small.
- **Keep the current name:** The name stays concise and the in-place modification is documented only in the function's docstring. No code changes elsewhere, but users who skim the name alone may be surprised when their input is mutated.
- **Keep the name but stop modifying the input:** The function gets a non-mutating semantics (copy the data internally). Cleanest API, but every call allocates a fresh copy of all localizations, so large datasets pay a time/memory cost on every run.

**Consumer input (@analysis, 2026-04-17):** Strong vote for non-mutating `cluster` with copy semantics. SMLMAnalysis's pipeline contract is "steps return new state"; mutation-by-default surprises users who pass the same SMLD to multiple steps or keep a pre-cluster reference. If performance pressure shows up later, add `cluster!` as a sibling — never mutate-and-return, that's the worst of both conventions. Emitters are lightweight structs, copy cost is negligible vs. the clustering itself.

**Technical detail:** `cluster(smld::BasicSMLD, cfg::AbstractClusterConfig) -> (smld_out, ClusterInfo)` currently mutates `emitter.id` on the input SMLD's emitters (emitters are mutable structs from SMLMData). The returned `smld_out` shares the input's emitter vector (or a filtered copy when `remove_unclustered=true`). Tests in `test_dbscan.jl:78-103` explicitly exercise and assert on this mutation. Downstream SMLMAnalysis uses this via `const DBSCANConfig = SMLMClustering.DBSCANConfig` and a thin `analyze(smld, cfg)` wrapper — renaming to `cluster!` requires updating that wrapper only. The docstring at `src/types.jl:60-78` documents the mutation; V2 in KB records the `(smld, ClusterInfo)` return convention (which the rename does not touch). If the answer is "stop modifying the input," the implementation cost is O(n) emitter-copy on every call; for a 10⁶-localization SMLD this is tens of MB and ~100 ms extra.

### Q4 — Should `ClusterInfo.cluster_sizes` be dataset-aware?

**From round:** 005

**Short question:** When the clustering is run on a multi-dataset file (imaging the same sample in several batches), the summary output reports one flat list of cluster sizes across the whole file. There is currently no way to tell from the summary which batch each cluster came from. Should the summary be changed to surface which batch each cluster belongs to?

**What each answer means:**
- **Yes — surface the batch per cluster:** The summary gains a second field pairing each cluster size with its batch number. Anyone reading the summary can immediately tell how clusters break down per batch. Any downstream code that already builds tables from the current summary needs to adapt to the added field.
- **No — keep the flat summary:** The summary stays simple. Users who need per-batch breakdowns recompute them from the per-localization labels on the output (which already carry the batch index).

**Consumer input (@analysis, 2026-04-17):** Flat `Vector{Int}` is fine for SMLMAnalysis's `step_summary` needs — it surfaces `n_clusters` and perhaps median size, nothing per-dataset. If a downstream user wants a per-batch breakdown they can group emitters by `(dataset, id)` themselves. Adding `Vector{Vector{Int}}` or `Dict` now is YAGNI; it's reversible if a real use shows up.

**Technical detail:** `ClusterInfo.cluster_sizes::Vector{Int}` currently concatenates local `1..K_local` sizes across datasets when `per_dataset=true`; `cluster_sizes[k]` is the size of the k-th cluster in visit-order across all datasets, and the mapping to `(dataset, id)` is lost at summary level. Tests at `test_dbscan.jl:115-144`, `test_hierarchical.jl:98-124`, and `test_voronoi.jl:97-129` implicitly assume this flat convention. Options if "Yes": (a) add `cluster_dataset::Vector{Int}` parallel to `cluster_sizes` (breaks positional constructor; tests use positional today, grep ClusterInfo usage in tests); (b) change to `Dict{Int, Vector{Int}}` (dataset → sizes) — more invasive; (c) switch to `Vector{NamedTuple{(:dataset, :id, :size), ...}}`. Option (a) is least invasive. ClusterInfo is re-exported by SMLMAnalysis so the change surfaces one layer up; docs field table in `src/types.jl:35-41` would update.

### Q5 — Is the Ward linkage + "cut in nanometers" combination on the Hierarchical backend a misleading default?

**From round:** 005

**Short question:** The hierarchical clustering backend asks users for a distance cut-off in nanometers, but its default linkage mode internally measures merges in a different unit (variance cost, not distance). That means passing "200 nm" under the default does not actually cut at 200 nm. Should the default be changed so the nanometer cut-off has its intended meaning, or should the field be renamed so it no longer implies nanometers?

**What each answer means:**
- **Change the default to a distance-based linkage (e.g. single-linkage):** The nanometer cut-off works as users expect straight out of the box. Advanced users who want variance-minimizing linkage opt in explicitly. SMLM-typical workflows gain a simpler mental model.
- **Rename the field so it no longer implies nanometers:** The cut-off name is generic (e.g. `cut_h`), and the docstring explains that the unit depends on the linkage. More honest, but every downstream caller has to learn what "cut height" means per linkage.
- **Split into two separate backend configs:** One for distance-based linkages (keeps `cut_nm`), one for variance-based linkages (different field name, different semantics). Cleanest but adds a type.
- **Leave as-is, document more loudly:** Keep the current defaults and field name, just hammer in the docstring that Ward's "nm" is a cost-unit, not a distance.

**Consumer input (@analysis, 2026-04-17):** Field name is misleading for Ward because the linkage "distance" isn't in nm — it's a variance-increase metric. Two clean options: (a) rename to `cut_threshold::Float64` with docstring calling out linkage-dependent units, or (b) add optional `n_clusters::Union{Int,Nothing}` alongside `cut_nm` so users can specify K directly (natural for Ward, awkward for single-linkage). @analysis would do both — rename for honesty, add `n_clusters` for ergonomics. If only one, prefer the rename. None of these change the SMLMAnalysis-consumed interface, so no blockers on that side.

**Technical detail:** `HierarchicalConfig` defaults to `linkage=:ward` and names its cut field `cut_nm::Float64` (converted to μm via `/1000.0` before passing to `Clustering.cutree(hc, h=cut_h)`). Ward linkage's dendrogram heights are sums-of-squares merging costs in μm² (not a distance), so `cut_nm=200.0` under Ward means "cut where the merging cost crosses 0.04 μm²," which users will not intuit. Single/complete/average linkages use Euclidean distance heights in μm where the "nm" label is correct. Current tests use `linkage=:single` throughout, so the mismatch is not exercised. KB V6 defines the `_nm` suffix convention for distance-valued fields. If the answer is "change default," only `src/backends/hierarchical.jl:40` changes; all tests remain passing (they all pass `linkage=:single` explicitly). If "rename," every test line that passes `cut_nm=` must update too.

---

## ANSWERED

<!-- Human has responded. Next round will process these and move them to PROCESSED. -->
<!-- Format keeps everything from OPEN and adds:
**Answer:** <human's direction>
**Answered:** YYYY-MM-DD
-->

(none)

---

## PROCESSED

<!-- Archived after a round applied the answer. Keep for audit trail. -->

### Q1 — Do we drop HDBSCAN from the package scope, or keep the slot open waiting for a lightweight library?

**From round:** 003

**Short question:** Should this clustering package ship without a hierarchical-density-based clustering backend, or should we keep that slot reserved for whenever a lightweight Julia implementation becomes available?

**What each answer means:**
- **Drop it:** The package ships with three backends (density-based, Voronoi, agglomerative hierarchical) and no hierarchical-density-based option. Anyone who wants the fourth algorithm wraps a Python one themselves.
- **Keep the slot open:** The package ships with three backends for now, but the fourth slot stays reserved in the roadmap and any future round that finds a suitable library picks it up.

**Technical detail:** D1 dead end in `KNOWLEDGE_BASE.md` records why HDBSCAN is blocked — `Clustering.jl 0.15.8` has no `hdbscan`, the only registered Julia HDBSCAN is in `HorseML.jl` (CUDA + NNlib + Zygote, too heavy), and `baggepinnen/HDBSCAN.jl` is GitHub-only and requires Python. The `@analysis` coordinator voted "drop it" to preserve the SMLMData-only dep footprint. If "drop it": remove Priority 3 from `STATUS.md` Future Priorities and note D1 as explicitly out-of-scope rather than BLOCKED. If "keep the slot open": leave Priority 3 as BLOCKED with a reference to D1; future rounds check periodically for new libraries.

**Answer:** Keep the slot open. Leave Priority 3 as BLOCKED referencing D1; future rounds can check periodically for new lightweight Julia HDBSCAN options.
**Answered:** 2026-04-17
**Processed:** Round 004 — STATUS.md Priority 3 annotated "Q1 resolved: keep the slot open and revisit when a suitable library appears." Priority 6 clarified that only HDBSCAN's tests remain pending (dependent on a future library).

### Q2 — How should the downstream analysis package bring this package in as a dependency?

**From round:** 003

**Short question:** Should the downstream analysis package pull this one in via a local development checkout, directly from its GitHub repository, or only after formal registry publication?

**What each answer means:**
- **Local development checkout:** Fastest iteration, both repos can evolve together; anyone building the downstream package needs the local checkout on their machine.
- **GitHub repository URL:** Slightly slower iteration but reproducible — commits are pinned by hash; works for anyone with internet access.
- **Wait for registry publication:** Cleanest for external users and reproducibility; delays downstream integration by however long registration takes.

**Technical detail:** `@analysis` (SMLMAnalysis) is waiting on this decision to wire `src/steps/cluster.jl` and add SMLMClustering to its `Project.toml`. Options map to: `Pkg.develop(path="/home/kalidke/julia_shared_dev/SMLMClustering")`, a `[sources]` entry with the GitHub URL, or a standard `Pkg.add("SMLMClustering")` after General-registry registration. Current state: SMLMClustering is version 1.0.0-DEV and not registered; lives at `/home/kalidke/julia_shared_dev/SMLMClustering`.

**Answer:** Hand the pick back to @analysis — SMLMAnalysis is the consumer, so whatever integration mode works best for them is the right one. Keith has already sent the handoff via agent-send; @analysis's reply will come into the clustering inbox. Next round should treat @analysis's response as the definitive answer and apply it (e.g. if they pick dev-path, make sure the package builds fine under `Pkg.develop`; if [sources], confirm the repo is pushed to GitHub; if wait-for-registration, queue a registration priority).
**Answered:** 2026-04-17
**Processed:** Round 004 — @analysis picked **option 2: `[sources]` with a GitHub URL, `rev="main"`**, mirroring the SMLMBaGoL precedent. Recorded as new Priority 7 in STATUS.md ("Push SMLMClustering to the JuliaSMLM GitHub org"). The actual `gh repo create` step needs Keith because it's an external action under the JuliaSMLM org; @analysis is waiting to be pinged with the URL + branch name.
