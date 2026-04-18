# Round 011 — Process ANSWERED Q3/Q4/Q5 (non-mutating cluster + cut_threshold/n_clusters)

**Date:** 2026-04-18
**Status:** done
**Priority worked on:** `QUESTIONS.md` ANSWERED items Q3 / Q4 / Q5 (processing takes precedence over Future Priorities per Phase 2 protocol). These were answered by Keith on 2026-04-17; Stop Condition was cleared same day.

## Hypothesis

Three separate API decisions converged in one round:

- **Q3 — non-mutating `cluster`.** Keith: "we don't modify our input; that is the interface convention used in SMLMAnalysis." @analysis's audit of the SMLMAnalysis step files confirmed universal non-mutation (one known exception on the `bagol.jl` diagnostic path, treated as bug-shaped). Implementation path: `deepcopy(smld.emitters)` at backend entry, rebind local `smld`, operate on the copy — per @analysis's option (b), which handles any `Emitter` subtype without per-type constructors.
- **Q4 — keep flat `cluster_sizes::Vector{Int}`.** No code change; mark PROCESSED.
- **Q5 — rename `cut_nm` → `cut_threshold` and add optional `n_clusters`.** The field name `cut_nm` was a lie under `:ward` (heights are variance costs, not distances). Rename honesty: `cut_threshold::Union{Float64,Nothing}=nothing` whose unit is linkage-dependent; `n_clusters::Union{Int,Nothing}=nothing` lets users sidestep the unit question for Ward. Exactly one must be set; both-set or neither-set → `ArgumentError`.

Going in: I expected the Q3 copy-semantics change to ripple into every test that asserted on `smld.emitters` post-call (DBSCAN, Hierarchical, Voronoi all had mirrored `labels written to emitter.id + remove_unclustered` test sets), and the Q5 rename to cascade to every `cut_nm =` call site. Unit tests already using `linkage = :single` would continue to validate at the same `cut_threshold` numeric values since nm→μm conversion is unchanged for distance-based linkages. Ward's pass-through (no division) is new code but exercised only through the new `n_clusters` path in this round's test additions.

## What was attempted

- **Phase 1 state read:** CLAUDE.md not present in project (global CLAUDE.md covers it), STATUS.md, the last two rounds (010 stopped + 009 idle-close), KNOWLEDGE_BASE.md (V1–V8, D1), QUESTIONS.md (three ANSWERED items), agent inbox.
- **Phase 1 inbox scan** (`~/.claude-agents/papers-64/inbox/clustering.jsonl`): two relevant messages on top of the backlog — @analysis 2026-04-17T13:53 confirming non-mutation is the SMLMAnalysis convention with the deepcopy-vs-rebuild implementation options, and @analysis 2026-04-18T09:29/09:30 ack'ing the P11/P12 scope additions (out of round's scope but shows coordinator is aligned and holding on wire-in until the GitHub push).
- **Q3 implementation.** Inserted a non-mutating guard at the top of every `cluster(smld, ::*Config)` method:
  ```julia
  smld = SMLMData.BasicSMLD(deepcopy(smld.emitters), smld.camera,
                            smld.n_frames, smld.n_datasets, smld.metadata)
  ```
  This rebinds the local `smld` to a fresh SMLD built over deep-copied emitters; all downstream writes to `emitter.id` and the `_build_output` call see only the copy. Three files touched: `src/backends/dbscan.jl`, `src/backends/hierarchical.jl`, `src/backends/voronoi.jl`. Tradeoff: one `deepcopy` per `cluster` call — O(n) emitter-copy, ~tens of MB for a 10⁶-localization SMLD (per Keith's Q3 technical detail estimate). Accepted as the correct cost of ecosystem-aligned semantics; if performance pressure appears, a `cluster!` variant can be added without changing the un-suffixed entry point (V9 scope).
- **Q3 docstring.** Rewrote the abstract-`cluster` docstring at `src/types.jl:66-80` to state explicitly "The input `smld` is **not modified**" and that each backend deep-copies before writing labels. Matches the V9 KB language exactly.
- **Q3 test updates.** For each backend's "labels written + remove_unclustered" testset I swapped:
  - Any assertion that ran over `smld.emitters` (the input) now runs over the output SMLD (`smld_keep.emitters`, `smld_out.emitters`, etc.) — where the labels actually live.
  - The pattern "rebuild a fresh SMLD because the previous call mutated ids in place" that DBSCAN/Hierarchical/Voronoi shared is gone: both calls reuse the same input `smld` since non-mutation makes that safe. Added `@test all(e -> e.id == 0, smld.emitters)` before/after each call as a regression anchor.
  - DBSCAN's `per_dataset label namespace` testset read `smld.emitters` to check per-dataset id ranges — switched to reading the returned `smld_out.emitters`. Same for Voronoi's equivalent testset. Hierarchical's got the same treatment plus a second call on the same `smld` (safe now).
- **Q4 implementation.** No code change required. PROCESSED marker only (Phase 4 file updates).
- **Q5 `HierarchicalConfig` rewrite.** Full rewrite of `src/backends/hierarchical.jl` to:
  - Swap `cut_nm::Float64` for `cut_threshold::Union{Float64,Nothing}=nothing`.
  - Add `n_clusters::Union{Int,Nothing}=nothing`.
  - Add mutual-exclusion validation: `xor(ct_set, nc_set)` or `ArgumentError`.
  - Branch the nm→μm conversion on linkage: `:ward` passes `cut_threshold` through unchanged (its unit is variance, roughly μm²); `:single`, `:complete`, `:average` divide by 1000 as before.
  - Dispatch `Clustering.cutree`: `cutree(hc, h=cut_h)` when `cut_threshold` drives, `cutree(hc, k=min(cfg.n_clusters, length(idxs)))` when `n_clusters` drives. The `min` clamp handles groups smaller than requested K (e.g. `n_clusters=10` on a 5-point group) cleanly without the underlying library throwing on out-of-range K.
  - Preserve the Q3 deepcopy block at the top.
- **Q5 test rewrite.** Full rewrite of `test/test_hierarchical.jl`:
  - Every `cut_nm =` → `cut_threshold =` (config constructor tests, all clustering testsets, argument-validation testset, `use_3d` error test, 3D path, empty SMLD, min_points filter).
  - New "n_clusters path (Ward linkage)" testset: three tight σ=10 nm blobs, `HierarchicalConfig(n_clusters=3, linkage=:ward, min_points=5, per_dataset=false)` → `n_clusters==3`, all 120 points clustered, zero noise. Demonstrates the Ward entry path with `n_clusters` rather than the broken-under-Ward `cut_threshold` numeric.
  - New error cases in "argument validation": `HierarchicalConfig(n_clusters=0)` → `ArgumentError`; `HierarchicalConfig()` (neither set) → `ArgumentError`; `HierarchicalConfig(cut_threshold=100.0, n_clusters=2)` (both set) → `ArgumentError`.
  - Non-mutation anchors on the "labels written" testset, as in DBSCAN and Voronoi.
- **Docs sync.** Updated `README.md` entry-point prose ("The call is **non-mutating**") and the `### Hierarchical` block (two code examples for the two paths, unit-convention paragraph rewritten). Updated `api_overview.md` cluster-entry-point "Side effects" block ("None. The input `smld` is not modified.") and the `HierarchicalConfig` field table, validation, unit caveat, and constructor examples.
- **KB updates.** Amended V6 with the linkage-dependent-unit caveat naming `cut_threshold` explicitly as the deliberate exception to the `_nm` suffix rule. Added V9 as the non-mutation invariant, with reference to Round 011 / Keith 2026-04-17 / @analysis's SMLMAnalysis-wide audit, and the deliberate asymmetry against the forthcoming passthrough `cluster_statistics` entry point.
- **Tests run:** `julia --project=. -e 'using Pkg; Pkg.test()'` → **175/175 passing in 29.1 s** (up from 159 in Round 008, +16 tests from the Q5 n_clusters path + mutual-exclusion errors + non-mutation anchors).

## What worked

- Non-mutating rebind at backend entry is surgical — three identical 4-line insertions, no helper touched, no invasive signature change. `_build_output`, `_group_by_dataset`, and the per-group label writes all operate on the local `smld` and don't know the copy happened.
- The "input unchanged after `cluster`" anchor tests caught one subtle thing the original code base had assumed: callers were building a fresh `_make_2d_smld(pts)` before every call to avoid "mutated ids." Removing those rebuilds in the updated tests simplifies each testset and is the most direct evidence the non-mutation claim holds — if it didn't, the second call on the same input would see labels from the first call and produce different cluster counts.
- Ward-path `n_clusters=3` test recovers exactly 3 clusters on three tight blobs with zero noise, confirming the `Clustering.cutree(hc, k=3)` branch works and the `min(cfg.n_clusters, length(idxs))` clamp is a no-op on normally-sized groups.
- Mutual-exclusion error tests pin the `xor(ct_set, nc_set)` behavior. Neither-set fires at the top of the method (no partial work done); both-set fires at the same check-point before the per-group loop.
- Per-linkage divide-by-1000 keeps all existing `:single`-linkage tests passing at identical `cut_threshold` numeric values (100.0, 200.0, 50.0) — the `:ward` branch is new code, but distance-based linkages are byte-identical to the pre-rename behavior.
- Docs + KB sync closes the interface-drift risk between source docstrings and top-level markdown. `api_overview.md` was the one most likely to rot silently; rewriting the `Side effects` and Hierarchical tables keeps Round 008's LLM-parseable claim valid.

## What failed

Nothing.

## Files changed

```
src/types.jl                            +3 -3    (cluster() docstring: non-mutating semantics)
src/backends/dbscan.jl                  +4 -0    (deepcopy rebind at entry)
src/backends/voronoi.jl                 +4 -0    (deepcopy rebind at entry)
src/backends/hierarchical.jl            full rewrite  (Q3 deepcopy + Q5 cut_threshold/n_clusters + docstring)
test/test_dbscan.jl                     +12 -14  (non-mutation anchors + switch to output smld)
test/test_voronoi.jl                    +10 -6   (non-mutation anchors + switch to output smld)
test/test_hierarchical.jl               full rewrite  (cut_nm→cut_threshold + n_clusters testset + mutual-exclusion errors + non-mutation anchors)
README.md                               +15 -7   (non-mutating prose + Hierarchical two-path examples + unit caveat)
api_overview.md                         +15 -9   (Side effects + Hierarchical field table + validation + caveat + constructor)
KNOWLEDGE_BASE.md                       +5 -1    (V6 caveat + V9 non-mutating invariant)
STATUS.md                               +N -M    (Current State refresh + Q3/Q4/Q5 PROCESSED callout + Round 011 row)
QUESTIONS.md                            +3 -0    (Q3/Q4/Q5 ANSWERED → PROCESSED)
rounds/round_011_process-q3-q5-answers.md   new
```

## Confidence

High. Full test suite 175/175 passing (29.1 s). Every assertion that previously read `smld.emitters` post-call either now reads the output SMLD (where the labels live) or explicitly checks the input was untouched (the non-mutation regression anchor). The Q5 rewrite preserves all pre-existing `:single`-linkage behavior byte-for-byte (same nm→μm conversion, same `cutree(hc, h=...)` call shape) while adding the Ward-safe `n_clusters` path with its own passing test. The KB V9 claim is the exact invariant the tests now assert on every backend.

## External consultations (if any)

None.

## Next steps

Remaining STATUS.md items split cleanly:

- **Priority 7 (GitHub org push)** — still blocked on Keith creating `JuliaSMLM/SMLMClustering.jl`. No round-worker action possible; @analysis is waiting on URL + branch to wire in via `[sources]`.
- **Priorities 11/12/13/14 (cluster_statistics + Hopkins)** — fully unblocked now that Q3/Q5 are landed and the interface is stable. P11 (scaffold the sibling interface) is the natural next pick and is concrete enough to be Sonnet-friendly once the pattern is established; the first P11 round itself has design content (new abstract type hierarchy, exports, passthrough-same-reference docstring) so Opus fits it better. P12+ (implement `HopkinsConfig`, then tests, then docs) is mechanical pattern-replication of the backend/test/docs triad we already have for clustering — Sonnet-appropriate rounds after P11 sets the interface.
- **Priority 3/6 (HDBSCAN)** — still D1-blocked; no new library landed since 2026-04-17.

Default: next round picks **Priority 11 — scaffold `cluster_statistics` interface** (abstract type + info struct + fallback + exports + KB V10). Keeping Opus on that round given the design content. No Sonnet opt-in written for next round.

## Questions posted

None. Q3 / Q4 / Q5 all moved ANSWERED → PROCESSED; no new OPEN items.

---

## Protocol notes (rare)

The previous round (010) was a `stopped` round triggered by the Stop Condition then in effect. That condition was cleared 2026-04-17 after Keith answered Q3/Q4/Q5 (STATUS.md Stop Conditions section was edited to `(none — ...)`). This round processed the ANSWERED items per the Phase 2 "Process answered questions" rule, which takes precedence over walking Future Priorities. The cadence table would normally have scheduled a review at Round 010; since that round was stopped (not counted as a productive round for cadence purposes), the next review now falls at Round 015 with the 2N KB audit at Round 020.
