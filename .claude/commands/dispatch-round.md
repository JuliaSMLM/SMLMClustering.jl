---
description: Spawn a fresh Claude session in a tmux window to run one /start-round with clean context. Prints a dashboard summary to this pane when the round completes.
---

# /dispatch-round — Fresh-Context Round Dispatcher

Thin wrapper around `.claude/scripts/dispatch-round.sh`. The script:

1. Reads config from the YAML frontmatter of `.claude/commands/start-round.md` (same project uses one config file, no separate TOML).
2. Checks the lock file — if held, skips. If held but the tmux window is gone (orphan), alerts.
3. Creates a new tmux window, launches a fresh Claude CLI session inside it, sends `/start-round`.
4. Polls the lock file in the background until it is released.
5. Kills the tmux window, reads state from disk, and prints a rich dashboard summary to this pane:
   - What the round just did (status, focus, key finding, commit)
   - Project health (test count, uncommitted changes, total rounds)
   - Upcoming priorities (next 3–5 items)
   - **Open questions awaiting your review** (from `QUESTIONS.md`)
   - Pointer to the next priority

## What to do

Run the script, read the dashboard it prints, then **post a synthesized summary to the user — not a raw dashboard dump.**

```bash
bash .claude/scripts/dispatch-round.sh
```

The dashboard is an audit/grep artifact for disk and `/loop` logs. The inline conversation deserves a human-readable digest. After the script exits, reply with:

- **2–3 sentences in your own words** on what the round did and why (synthesized from the "What this round did" block and the commit message — do not copy verbatim).
- **Any OPEN items** from the "Open questions" block, transcribed as-is. Say "(none)" if empty.
- **Next 3–5 priorities**, copied from the "Upcoming priorities" block with the status tags intact.
- The round's **status** (done / partial / aborted) and **commit SHA** as a one-line footer.

Skip the double-line separators, the "Project health" metrics, and the sectional headers — those are dashboard chrome, not signal. Do not do round work yourself.

## After posting the digest — refresh dashboard views (if enabled)

If `.claude/dashboard/` exists in the project root, the project has the `round-dashboard` skill active. After posting your digest reply, follow the instructions in `.claude/dashboard/regenerate.md` to refresh:

- `.claude/dashboard/status-view.md` — compact status summary for the top-left tmux pane
- `.claude/dashboard/questions-view.md` — numbered Q1/Q2/Q3 list for the top-right tmux pane

The regenerate prompt specifies the exact output shape. Cost: two short `Edit` calls. The tmux watchers in the dashboard panes pick up the changes within ~1 second of your edits finishing.

If `.claude/dashboard/` is absent, skip this step — the project is running plain round-init without the dashboard layer. No-op, no behavior change.

## Usage

Manually: `/dispatch-round`

With a loop: `/loop 15m /dispatch-round` — fires every 15 minutes. If the previous round is still running (lock held), the new fire is skipped, so there's no overlap. The dispatch dashboard prints between rounds so you always see what just happened, what's next, and what needs your attention.

## Answering questions (dashboard mode only)

If `.claude/dashboard/` exists and the user's message starts with `q<N>` (case-insensitive: `q1`, `Q2`, `q3 — yes, do X`), interpret it as an answer to the Nth question in the current `.claude/dashboard/questions-view.md`. Procedure is defined in `.claude/dashboard/regenerate.md` under "Answer shorthand". Summary:

1. Read `.claude/dashboard/questions-view.md` and match `q<N>` to a `## Q<N>.` heading.
2. Edit the full `QUESTIONS.md` at repo root: move that item from `## OPEN` to `## ANSWERED`, appending the user's answer and today's date.
3. Regenerate `.claude/dashboard/questions-view.md` (the answered question disappears).
4. Confirm in one sentence: `"Q<N> moved to ANSWERED. Next round will process."`

## Why this exists

Running a dev prompt in one long session degrades: after ~15–20 rounds, accumulated context pollutes judgment. Each round in its own fresh tmux session with clean context avoids that. The state files (`STATUS.md`, `KNOWLEDGE_BASE.md`, `QUESTIONS.md`, `rounds/`) are the institutional memory that makes fresh-context rounds viable — nothing carries over in context, everything persists on disk.
