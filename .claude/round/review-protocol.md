# Review Protocol + KB Audit

Structured handling of `/review-code` output and the every-2N-review-round KB audit. Both follow fixed procedures so the model does not improvise.

---

## Part 1 — Code Review (every `review_cycle_rounds`)

When Phase 2 determines this is a review round, run `/review-code` and process its output with the following protocol. Do not silently drop any finding.

### For each review finding

Classify into one of four buckets:

| Bucket | Meaning | Action |
|---|---|---|
| AGREE — trivial | You agree and it's a small fix (typo, unused import, minor refactor) | Apply immediately in this round |
| AGREE — substantial | You agree but it's larger than a trivial fix (architecture change, broad refactor) | Add as a new item to `Future Priorities` in `STATUS.md` with reference to the review round |
| DISAGREE | You think the finding is wrong, or premature, or based on a misreading | Post to `QUESTIONS.md` as an `OPEN` item with your reasoning. Do NOT just ignore. |
| DEFER | Agree it's worth thinking about but needs human input | Post to `QUESTIONS.md` as `OPEN` |

### Review round output

Record in the round file (`rounds/round_NNN_review-at-NNN.md`):

```
## Review findings

- [AGREE-trivial] <finding> — fixed in <commit hash>
- [AGREE-substantial] <finding> — added as Priority <N> in STATUS.md
- [DISAGREE] <finding> — see QUESTIONS.md Q<n>. Reasoning: <...>
- [DEFER] <finding> — see QUESTIONS.md Q<n>
```

The review round itself replaces normal priority work for this round.

---

## Part 2 — KB Audit (every 2 × `review_cycle_rounds`)

Every second review round is also a KB audit round. With default `review_cycle_rounds: 5`, this fires at round 10, 20, 30, ...

The audit is **mechanical, not exploratory**. One pass through the KB with a fixed decision taxonomy. It should take a few minutes, not half an hour.

### Inputs

- `KNOWLEDGE_BASE.md` — all V-entries, D-entries, prior A-entries
- Round files since the last audit (up to 10 rounds)
- The `/review-code` output from this round (preceding the audit in the same round)

### Procedure

**Step A — Walk each V-entry and D-entry once.** Classify into exactly one outcome:

| Outcome | Meaning | Action |
|---|---|---|
| **Confirm** | Still true; evidence from rounds since last audit supports it | Leave as-is |
| **Contradict** | Recent findings directly disprove or invalidate it | Mark `RETIRED` inline. Add back-reference to the round file with the contradicting evidence. DO NOT DELETE. |
| **Narrow** | Still partially true, but the original claim was too broad | Edit the entry to the narrower scope. Note the narrowing with a back-reference. |
| **Supersede** | Replaced by a newer, better entry (often a D that became a V) | Mark `SUPERSEDED by V<n>` (or `D<n>`). Leave original in place. |

**Step B — Consolidate.** Walk the list again. If two entries say substantively the same thing, merge into one and mark the duplicate `MERGED into V<n>` with a back-reference.

**Step C — Emit a new A-entry.** Append to `KNOWLEDGE_BASE.md` under the Audits section:

```markdown
### A<next> — Audit at Round NNN

Reviewed: V1–V<max>, D1–D<max>
Confirmed: <comma-separated list>
Retired: <list with back-references to contradicting round files>
Narrowed: <list with before/after scope summary>
Superseded: <list>
Merged: <list>

Summary: <X of Y entries confirmed unchanged>. <Key takeaway, one sentence>.
```

**Step D — Record in the round file.** The round file for an audit round mentions the A-entry id and lists which entries moved status (same content as in the A-entry, but in the round file for round-level auditability).

### What the audit is NOT

- NOT a re-validation from scratch of every decision
- NOT a re-read of every round file in detail
- NOT an exploration of "are there better ways" for each entry
- NOT a ritual — if the KB is small or stable, most entries will confirm unchanged and the audit finishes quickly

### Why this cadence

After ~30–40 rounds, stale V-entries start getting read in Phase 1 every round, consuming token budget and occasionally leading the worker down paths that were already known to be wrong. The audit prevents compounding staleness without costing much per cycle.
