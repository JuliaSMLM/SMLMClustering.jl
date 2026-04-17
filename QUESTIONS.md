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

(none)

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
