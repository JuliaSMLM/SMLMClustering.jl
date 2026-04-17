# Dashboard View Regeneration

This file is read by the parent agent after each `/dispatch-round` completes. It governs the two short summary view files displayed in the tmux dashboard panes. This file is project-editable — tune the formatting if you want different emphasis.

The parent agent should follow this procedure once per round, **after** posting its digest reply to the user in the conversation. Cost: two short `Edit` tool calls. Negligible overhead.

---

## Step 1 — Regenerate `.claude/dashboard/status-view.md`

Read `STATUS.md` at the project root. Write `.claude/dashboard/status-view.md` with exactly this shape (no extra sections, no header changes):

```markdown
# {{PROJECT_NAME}}

## Project Overview

<2–3 sentences. What the project is scientifically — the research question and the shape of the deliverable. Static-ish; changes only when the project's fundamental scope changes, not round-to-round. Plain-English, avoid project-internal jargon like "Priority N" or "Q<n>" — name the work semantically. Derived from `CLAUDE.md` and the current `## Current State` body of `STATUS.md`.>

## Development Overview

<2–3 sentences. What phase the development effort is in right now: headline result so far, what is left, what is blocking. Rotates per round. Same plain-English rule as Project Overview — replace "Priority 4" with the semantic name of that priority, replace "Q3" with the semantic name of that question.>

## Task List

### Recent

| N | Focus | Status |
|---|-------|--------|
| <NNN> | <focus> | <done/partial/aborted> |
| <NNN> | <focus> | <done/partial/aborted> |
| <NNN> | <focus> | <done/partial/aborted> |

### Last round

<2–3 sentences. Plain-English first sentence: what was accomplished and why it matters — no file paths, no commit SHAs, no tool names, no priority numbers. Then a single `**Technical detail:**` paragraph with the specifics (file paths, commit SHA, tool or function names, key numbers). Rule: if you removed the `**Technical detail:**` paragraph, the opening should still answer "what did this round do?" — if it doesn't, the framing is wrong.>

<!-- ACTIVE -->

### Next

<Top N highest-severity open priorities. Walk `Future Priorities` top to bottom. Skip items marked `DONE`, `COMPLETED`, `BLOCKED`, `~~strikethrough~~`, `✅`, or `EFFECTIVELY DELIVERED`. Among the remaining, find the highest severity tag present (`CRITICAL` > `HIGH` > `MEDIUM` > `LOW`) and list up to 3 at that level. If fewer than 3 exist at the highest level, fill from the next level down. If 5 items all tie at `MEDIUM`, show all 5 — the user wants to see what is urgent, not what is next in file order. If no explicit severity tags are present, show the first 3 non-skipped items.>

1. **[{severity}]** <priority text> — <status>
2. ...
```

Filter rules for "Recent":
- Take the last 3 rows of the Round History table. If fewer than 3 rounds have run, show what exists and pad with `| — | — | — |` to keep the table well-formed.

**Preserve the `<!-- ACTIVE -->` marker line.** The dashboard's render-view splice-mode watcher replaces that marker with the contents of `/tmp/<slug>-active.md` when a round is running (showing the active pick between `### Last round` and `### Next`) and drops it when no round is active. The marker is an HTML comment so it's invisible to mdcat when no splice happens. Keep it on its own line, exactly as shown above.

**No-jargon check.** Before saving `status-view.md`, scan the `## Project Overview` and `## Development Overview` paragraphs and the `### Last round` opening sentence. If any of them contain `Priority N`, `Priority #`, `Q<digit>`, or `P<digit>` (the internal shorthand), rewrite that phrase with the semantic name of the work. The view file is the user's primary at-a-glance surface; it must read like prose, not like a ticket queue.

---

## Step 2 — Regenerate `.claude/dashboard/questions-view.md`

Read `QUESTIONS.md` at the project root. Write `.claude/dashboard/questions-view.md` with this shape:

```markdown
# Open Questions

## Q1. <short question ending in ?>

**What each answer means:**
- **Yes / Option A:** <outcome>
- **No / Option B:** <outcome>

## Q2. <short question ending in ?>

**What each answer means:**
- **Yes / Option A:** <outcome>
- **No / Option B:** <outcome>

...
```

Rules:
- Renumber from `Q1`, `Q2`, `Q3` regardless of how the OPEN items are anchored in `QUESTIONS.md`. The display numbering is for the user's reference; the source file's anchors don't change.
- Drop `**From round:**` and `**Technical detail:**` blocks — those are for workers, not the human. Keep `**Short question:**` content as the heading line, and `**What each answer means:**` verbatim.
- If the OPEN section is empty, write:

```markdown
# Open Questions

(none)
```

---

## Answer shorthand — `q<N>`

When the user replies in the parent pane with a message starting `q<N>` (case-insensitive: `q1`, `Q2`, `q3` …), interpret it as an answer to the Nth question in the current `questions-view.md`:

1. Read `questions-view.md` to identify which `## Q<N>.` heading the user means.
2. Match it against the OPEN items in the full `QUESTIONS.md` at repo root (titles should match; if ambiguous, ask the user to disambiguate).
3. Edit `QUESTIONS.md`: move the item from `## OPEN` to `## ANSWERED`, appending:
   ```markdown
   **Answer:** <user's answer text>
   **Answered:** <today's date YYYY-MM-DD>
   ```
4. Regenerate `questions-view.md` (Step 2 above) so Q<N> disappears from the display.
5. Confirm to the user in one sentence: `"Q<N> moved to ANSWERED. Next round will process."`

The next `/dispatch-round` will PROCESS the ANSWERED item per the round-init protocol (Phase 2) and update STATUS.md / KB accordingly.
