---
description: Orchestrator-side dispatcher for the fresh-context round system. Pre-flight queue check (with auto-stop), spawns the round worker as a forked Agent, refreshes the dashboard. No shell script — the orchestrator drives everything via Agent tool calls.
---

# /dispatch-round — Orchestrator-Side Round Dispatcher

You (the orchestrator) are the long-running session that runs `/loop /dispatch-round`. Your accumulated chat with the developer + project narrative is **context** (signal). Tool work — file reads, edits, commits, dashboard refresh — is **pollution** that you must keep out of your context by delegating to forked sub-agents.

This file is what you do each time `/loop` fires `/dispatch-round`.

**Mental model:** you orchestrate; forks execute. The chat survives; the tool noise doesn't.

---

## 1. Pre-flight (small fork — queue check + model recommendation)

Spawn a forked Agent to read state and recommend the round shape. This keeps STATUS.md / QUESTIONS.md text and frontmatter parsing out of your own context.

```
Agent({
  description: "Round NNN pre-flight",
  prompt: "Read .claude/commands/start-round.md (frontmatter only — note model_mode and lock_slug),
           STATUS.md (Future Priorities + Stop Conditions sections),
           QUESTIONS.md (OPEN section). Then return a compact summary in this exact form:

           NNN: <next zero-padded round number from ls rounds/round_*.md sort tail -1, +1>
           ACTIONABLE: <count of priorities under Future Priorities that are TODO or IN PROGRESS,
                       not BLOCKED, not DONE/COMPLETED, not strikethrough/✅,
                       AND have all three of add:, test:, doc: lines>
           UNDER_SPECIFIED: <count of priorities that would be actionable but lack the trio>
           OPEN: <count of OPEN questions in QUESTIONS.md>
           STOP: <text of any active Stop Condition, or 'none'>
           TOP: <text of first actionable priority, or '(none)'>
           NEEDS_OPUS: <true if TOP carries the [needs_opus] tag, else false>
           MODEL_REC: <opus|sonnet — apply model_mode rule:
                       opus_always → opus
                       sonnet_always → opus if NEEDS_OPUS else sonnet
                       bias_opus → sonnet ONLY if TOP is clearly mechanical
                                  (pattern replication, test/doc additions to existing code,
                                   docstring/README pass, rename/deprecation cleanup, no semantic change);
                                  otherwise opus
                       bias_sonnet → opus if NEEDS_OPUS or TOP is clearly judgment-heavy
                                    (ambiguous criteria, cross-module reasoning, design decision
                                     framing, prior round was partial/aborted on this priority,
                                     priority text asks a question);
                                    otherwise sonnet>

           Do NOT modify any file. Do NOT execute round work. Read-only."
})
```

The fork returns the summary; you read it and proceed.

---

## 2. Auto-stop check

**If `STOP` is non-`none`:** the project is paused.
- Post a one-liner to the developer: `Round paused — Stop Condition: <text>. /loop continuing to check.`
- Skip the round (do NOT spawn the worker).
- Schedule the next /loop tick if dynamic mode (cron mode auto-fires).
- End this dispatch.

**If `ACTIONABLE == 0` AND `OPEN == 0`:** the queue is empty and there's nothing for a human to weigh in on. Terminate the loop:

1. Post to the developer:
   ```
   Queue empty — /loop terminating. No actionable priorities, no open questions.
   Add a feature-trio priority to STATUS.md (add: / test: / doc:) and re-arm with
   `/loop 15m /dispatch-round` when ready.
   ```
2. **Cron mode:** call `CronList`, find your own `/dispatch-round` cron entry, call `CronDelete` on it.
3. **Dynamic mode (`/loop` with no interval):** simply do not call `ScheduleWakeup` for this dispatch.
4. Spawn one final dashboard-refresh fork (Step 5 below) so the panes show queue-empty state.
5. End. Do not proceed to Step 3.

**If `ACTIONABLE == 0` AND `OPEN > 0` AND `UNDER_SPECIFIED == 0`:** waiting for human input. Skip the round but keep /loop alive — the next ANSWERED question will become the next round's seed.
- Post: `Round skipped — waiting on <count> OPEN question(s). /loop continuing.`
- Schedule next tick if dynamic mode.

**If `ACTIONABLE == 0` AND `UNDER_SPECIFIED > 0`:** there are priorities but they're missing the feature trio (`add:`/`test:`/`doc:`). The round worker would skip them too. Spawn a small fork to post OPEN questions asking the developer to specify each:
```
Agent({
  description: "Surface under-specified priorities",
  prompt: "Read STATUS.md Future Priorities. For each TODO/IN PROGRESS item missing
           one of add:/test:/doc: lines, post an OPEN question to QUESTIONS.md with
           the four-block shape (Short question ends in '?', What-each-answer-means,
           Technical detail). Skip if an identical question is already OPEN.
           Return: count of new OPEN questions added.",
})
```
Then post: `<count> priorities need add/test/doc spec — see QUESTIONS.md. /loop continuing.` Skip this round.

**Otherwise:** proceed to Step 3.

---

## 3. Round fork (the actual work)

Spawn the round worker. Pass the round number, the chosen model, and a short directive — the fork inherits your accumulated chat, so the project narrative is already there.

```
Agent({
  description: "Round NNN — <truncate TOP to ≤60 chars>",
  model: <MODEL_REC from pre-flight: 'opus' or 'sonnet'>,
  prompt: "Run a fresh round per .claude/commands/start-round.md (Phases 0–5).

           Round number: NNN
           You are running as: <opus|sonnet>

           Pick the priority via Phase 2 — the pre-flight identified TOP as
           '<TOP text>' but defer to your own Phase 2 walk if ANSWERED questions
           or stop conditions changed since pre-flight.

           Append progress to /tmp/<lock_slug>-worker-progress.md as you go
           (the user's dashboard worker pane is tailing it).

           Return the round summary in the format specified by Phase 5."
})
```

The Agent call blocks until the fork returns, but **you can chat with the developer in the meantime** — fork tool output never enters your context. Type-while-fork-runs is the canonical queue-while-running surface.

While the fork runs, expect inbound messages from the developer:

| Message shape | Action |
|---|---|
| `q<N> <answer>` (e.g. `q2 yes use Reactant`) | Spawn small fork to process the answer (move OPEN→ANSWERED in QUESTIONS.md, regenerate questions-view.md). See "Answer shorthand" below. |
| `add priority X to STATUS` (or similar) | Spawn small fork: read STATUS.md, append a feature-trio entry (or post OPEN if the developer's text under-specifies). |
| `abort the current round` | Send an interrupt to the in-flight Agent fork. Round fork should reach Phase 4 and close as `aborted`. |
| Any other discussion | Reply inline. Don't fork. The conversation is what the orchestrator's context is *for*. |

If the developer asks for project facts that require reading files (e.g. "what does priority 3 say verbatim?"), spawn a small fork rather than reading the file yourself. Trivial single-line edits to state files are fine inline; **multi-file reads always fork.**

---

## 4. On round-fork return

Synthesize the digest **in your own words** and post it to the developer. Do not dump the fork's return verbatim.

Format:

- **2–3 sentences** on what the round did and why (what changed, why it matters — synthesized from the fork's Status / Focus / Key finding lines plus your own project context). Avoid jargon and internal shorthand ("Priority 4", "Q3") — name the work semantically.
- **Any OPEN items** added this round, transcribed from the fork's "OPEN questions added" line. Say `(none)` if empty.
- **Status + commit footer**: one line — `<done|partial|aborted> · <commit SHA>`.

Skip raw round-file contents, raw STATUS.md sections, and project-health metrics. Those are on disk; the dashboard refresh below surfaces them.

---

## 5. Refresh the dashboard view files

If `.claude/dashboard/` exists, the project has the round-dashboard layer active. Spawn a small fork to regenerate the view files:

```
Agent({
  description: "Refresh dashboard views",
  prompt: "Read .claude/dashboard/regenerate.md and follow its instructions exactly.
           Regenerate .claude/dashboard/status-view.md from STATUS.md and
           .claude/dashboard/questions-view.md from QUESTIONS.md per the shapes
           in regenerate.md. Two Edit tool calls expected. Return: 'done'."
})
```

If `.claude/dashboard/` is absent, skip — the project is running plain round-init without the dashboard. No-op.

The tmux watchers in the dashboard panes pick up the changes within ~1 second.

---

## 6. Schedule the next tick (dynamic mode only)

If you're running `/loop` in dynamic mode (no fixed interval), call `ScheduleWakeup` to fire `/dispatch-round` again. Default interval: 15 minutes. The auto-stop logic in Step 2 will terminate the loop when the queue is genuinely dry.

If you're running `/loop` in cron mode (`/loop 15m /dispatch-round`), the next tick fires automatically — do nothing.

---

## Answer shorthand — `q<N>` references

When the developer replies in your pane with a message starting `q<N>` (case-insensitive: `q1`, `Q2`, `q3 — yes, do X`), interpret `q<N>` as an answer to the Nth question in the current `.claude/dashboard/questions-view.md`.

Spawn a small fork:

```
Agent({
  description: "Process q<N> answer",
  prompt: "Read .claude/dashboard/questions-view.md to identify which OPEN question
           Q<N> refers to. Match it against QUESTIONS.md OPEN items by title.
           Edit QUESTIONS.md: move the item from ## OPEN to ## ANSWERED, appending
           '**Answer:** <user's text>' and '**Answered:** YYYY-MM-DD'. Then
           regenerate .claude/dashboard/questions-view.md per
           .claude/dashboard/regenerate.md. Return: the title of the Q that moved."
})
```

When the fork returns, post a one-line confirmation to the developer: `Q<N> moved to ANSWERED. Next round will process.`

If `.claude/dashboard/` is absent, fall back to matching `q<N>` against the Nth OPEN item in `QUESTIONS.md` directly.

---

## Why the orchestrator never reads files itself

The orchestrator session can run for weeks. Each `Read` you do directly bakes that file's content into your context for the rest of the session — even after the file changes on disk. That's pollution.

Forks invert the cost: tool noise lives in the fork's ephemeral context, returns to you as a short summary, and the next time you need that information you fork again with the latest disk state. The orchestrator's context stays the developer-discussion narrative, which is what makes future rounds richer (the round fork inherits it).

The exception is **trivial inline edits** to state files where you already know the change (e.g. the developer says "add priority X with add: foo, test: bar, doc: baz" — that's one Edit call to STATUS.md, no reads required, fine to do inline). When in doubt: fork.
