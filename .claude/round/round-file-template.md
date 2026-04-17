# Round File Template

Each round writes `rounds/round_NNN_<slug>.md` using this structure. `NNN` is zero-padded (001, 002, ...). Slug is a short kebab-case description of the focus.

---

```markdown
# Round NNN — <focus>

**Date:** YYYY-MM-DD
**Status:** done | partial | aborted | stopped
**Priority worked on:** <reference to priority in STATUS.md, e.g. "Priority 4: Remove train_zygote!">

## Hypothesis

One paragraph. What did we expect to be true going into this round? What outcome did we think the work would produce?

## What was attempted

Bulleted list. For each approach or sub-step:
- What was tried
- What files/functions touched
- Result (worked / failed / partial)

## What worked

The specific things that now function. Concrete, with test references or behavioral evidence. If nothing worked, write "Nothing — this was an abort/stop round. See 'Why it stopped' below."

## What failed

The specific things attempted that did not work. For each:
- What was tried
- How it failed (error message, unexpected behavior)
- Whether the failure changed our understanding

## Files changed

```
src/path/to/file.jl          +23 -5
test/test_whatever.jl        +12
STATUS.md                    +3 -1
KNOWLEDGE_BASE.md            +8
```

## Confidence

High | Medium | Low — on the claim that this round advanced the stated priority correctly.

Brief justification (one sentence). "Full test suite passes and the previously-failing integration test now covers the new code path" is high confidence. "Targeted test passes but I haven't run the full suite" is medium. "Change compiles but I couldn't verify behavior because X" is low.

## External consultations (if any)

If `/ask-codex` or `/second-opinion` was used, record each:
- Source: codex | second-opinion
- Question asked
- Outcome applied (or not applied, with reason)

## Next steps

What should the next round pick up? Can be:
- "Continue priority N, next sub-step is X"
- "Priority N now DONE, next unblocked priority is M"
- "Priority N now BLOCKED, see QUESTIONS.md item Q<n> for human judgment"
- "Aborted — 3-strike dead end recorded as D<n>. Priority re-prioritized in STATUS.md."

## Questions posted

If this round added to `QUESTIONS.md`, list the Q-ids and one-line summaries here for at-a-glance reference.

---

## Protocol notes (rare)

Only include this section if the round did something structural: KB audit results, stop condition triggered, human override applied from QUESTIONS.md, etc. Otherwise omit.
```
