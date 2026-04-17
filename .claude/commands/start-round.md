---
# --- Dispatcher config (auto-filled by /round-init; edit if paths move) ---
project_name: "SMLMClustering"
project_dir: "/home/kalidke/julia_shared_dev/SMLMClustering"
lock_path: "/tmp/smlmclustering-round.lock"
window_name: "smlmclustering-round"
launcher: "ccdspy"
startup_delay_s: 10

# --- Protocol policy (three knobs — see reference/rationale.md) ---
scope_mode: run_until_done          # iteration_cap | run_until_done
lock_strategy: manual_reclamation   # pid_based | manual_reclamation
review_cycle_rounds: 5              # KB audit runs every 2nd review round (forced)
iteration_cap: 8                    # used only if scope_mode == iteration_cap
---

# /start-round — Fresh-Context Round Worker

Runs inside a freshly-spawned Claude session. Does **one** atomic unit of work: reads state from disk, picks a priority, executes it, commits, updates state files, releases the lock.

All state persists via disk — nothing carries over between rounds except what is written to `STATUS.md`, `KNOWLEDGE_BASE.md`, `QUESTIONS.md`, and `rounds/round_NNN_*.md`.

**Context-sizing principle.** Each round is a fresh session with a 1M context window and an in-round prompt cache. Load relevant state eagerly — the cache makes re-references inside a round essentially free, and rationing context to "stay small" just forces reactive disk re-reads in Phase 3. Productive rounds typically land at 200k–400k of loaded context; do not self-police below that. Use the ceiling only if the priority genuinely requires it.

---

## Phase 0 — Lock

**Acquire the round lock atomically.** If another round is in progress, exit immediately.

```bash
LOCK=$(grep '^lock_path:' .claude/commands/start-round.md | awk '{print $2}' | tr -d '"')
if ! (set -o noclobber; echo "$$" > "$LOCK") 2>/dev/null; then
  echo "ROUND SKIPPED — lock held at $LOCK"
  exit 0
fi
trap 'rm -f "$LOCK"' EXIT
```

**Lock strategy** (from frontmatter `lock_strategy`):

- `pid_based` — write PID to lock file. If lock exists but PID is dead (`kill -0` fails), reclaim. Use for short deterministic rounds.
- `manual_reclamation` — if lock exists, exit immediately with no reclamation. Rounds may legitimately run for hours (training jobs, long fits). Staleness is the dispatcher's problem (window-gone-but-lock-exists → orphan cleanup) or the human's problem (explicit `rm`). This is the default.

**On any exit (normal, error, abort) the lock MUST be released.** Phase 4 writes; the trap is the backstop.

<!-- PROJECT-SPECIFIC: phase-0-extras — preserved across rescaffolds -->
<!-- Add project-specific lock behavior here, e.g. clearing /tmp/*.git.lock -->
<!-- /PROJECT-SPECIFIC -->

---

## Phase 1 — State Read

**Strict order. Read generously — bounded by priority-selection needs, not context cost.** The fixed prefix (state files in this order) primes the in-round cache; rationing here just shifts cost into Phase 3 as reactive re-reads.

1. `CLAUDE.md` (project instructions — this is the de-facto immutable anchor)
2. `STATUS.md` — full read (Current State, Active Threads, Future Priorities, Round History, Stop Conditions)
3. Latest 1–2 files in `rounds/` (most recent round files only)
4. `KNOWLEDGE_BASE.md`:
   - **Validated Decisions**: read in full
   - **Dead Ends**: read in full (you want to recognize a resembling priority before starting work, not mid-Phase-3)
   - **Audits**: read the most recent two
5. `QUESTIONS.md` — full read (check for `ANSWERED` items to process this round)

**Stop reading after this.** Additional file reads belong to Phase 3, driven by the chosen priority — not exploratory browsing.

<!-- PROJECT-SPECIFIC: phase-1-extras -->
6. **Agent inbox read.** Pull any cross-agent feedback that arrived since the last round. @analysis (SMLMAnalysis) coordinates this package's interface and may flag drift or request changes between rounds.

   ```bash
   SESSION=$(tmux display-message -t "$TMUX_PANE" -p '#{session_name}')
   AGENT=$(tmux display-message -t "$TMUX_PANE" -p '#{window_name}')
   INBOX="$HOME/.claude-agents/$SESSION/inbox/$AGENT.jsonl"
   [[ -f "$INBOX" ]] && tail -50 "$INBOX"
   ```

   Read every message — do not truncate based on "looks stale." A priority-shift request from @analysis must weight Phase 2 selection even if it arrived several rounds ago. If a message raises a design question, convert it into an OPEN item in `QUESTIONS.md` during Phase 4; do not answer on @analysis's behalf.
<!-- /PROJECT-SPECIFIC -->

---

## Phase 2 — Priority Selection

**Check stop conditions first.** Read the `Stop Conditions` section of `STATUS.md`. If any condition is met, halt with `ROUND STOPPED: <condition>` and proceed directly to Phase 4 (close with aborted-round round file).

**Process answered questions.** For each `QUESTIONS.md` item in `ANSWERED` status: apply the human's guidance, update `STATUS.md` / `KNOWLEDGE_BASE.md` as directed, mark the question `PROCESSED`.

**Walk Future Priorities top-to-bottom.** Pick the first item that is:
- `TODO` or `IN PROGRESS` (not `BLOCKED`, not `DONE`)
- `HIGH` or `CRITICAL` if severity tags are used
- Not blocked by an active thread
- Not resembling a KB Dead End (check Dead End titles)

**Cadence rules** (forced defaults):

| Round number | Extra action |
|---|---|
| Every Nth round where N = `review_cycle_rounds` | Run `/review-code`, process findings via `review-protocol.md` (AGREE/DISAGREE per finding) |
| Every 2N round (i.e. every second review round) | Run `/review-code` AND the KB audit procedure in `.claude/round/review-protocol.md` |

A review or audit **replaces** normal priority work for this round. Record as focus `Review at Round NNN` or `Review + KB Audit at Round NNN`.

**Announce the pick to the dashboard.** Immediately after priority selection (before starting Phase 3), write a small markdown file naming what you just picked. The dashboard's status pane splices this file into the status view so the active pick shows live, independent of round-close state.

```bash
printf '### Active pick\n\nRound NNN · Priority <id> · <short semantic description>\n' > /tmp/smlmclustering-active.md
```

Examples (the line after the header):

- `Round 011 · Priority 4 · write-up (Mann-Whitney tables + methods text)`
- `Round 012 · Priority 6 · figure set (per-condition distributions panel)`
- `Round 013 · Review at Round 015 · /review-code + classify findings`

Rules:

- Single `Bash` tool call with `printf '...' > /tmp/smlmclustering-active.md`. **Do NOT use the `Write` tool.** The harness classifies `.claude/` and some adjacent paths as sensitive; the Write-tool path prompts even under bypass-permissions. `/tmp/` via Bash `printf` / `echo` goes through the Bash tool's permission pathway (already approved for round workers) and does not prompt.
- Exact shape: `### Active pick` header, blank line, one content line ≤100 chars (no further markdown decoration on the content line).
- Phase 4 clears the file; see the close sequence.
- Write this even if `.claude/dashboard/` does not exist on the project. The file is cheap, keeps the protocol uniform across round-init projects, and costs nothing when unused.

<!-- PROJECT-SPECIFIC: phase-2-extras -->
<!-- Project-specific priority-selection overrides (rare). -->
<!-- /PROJECT-SPECIFIC -->

---

## Phase 3 — Execute

**Do the work.** Read `.claude/round/anti-patterns.md` and `.claude/round/project-anti-patterns.md` if you need a refresher on what not to do. Read `.claude/round/failure-modes.md` if you hit an unexpected failure.

**Hard rules** (forced, no override):

1. **Max 3 failed attempts** on the same obstacle. On the 3rd failure: stop, document the obstacle as a new Dead End in `KNOWLEDGE_BASE.md`, move on (or end the round).
2. **No scope creep.** If you discover something worth doing that isn't the current priority, add it to `Future Priorities` in `STATUS.md` — do not act on it this round.
3. **External consultation budget: 2 per round.** `/ask-codex` and `/second-opinion` are allowed only AFTER an in-tree attempt has failed. Each consultation must be recorded in the round file with the consulted source.
4. **Never modify** `CLAUDE.md`, `.claude/commands/start-round.md`, or canonical files in `.claude/round/` (anti-patterns, failure-modes, round-file-template, review-protocol). Editable files: `STATUS.md`, `KNOWLEDGE_BASE.md`, `QUESTIONS.md`, new round files, project source, `.claude/round/project-*.md`.

**Scope bounding** (from `scope_mode`):

- `run_until_done` — work until the priority is complete or you cannot make headway. Do not stop because the round feels long. Fresh-context sessions keep the relevant files in KV cache; ending early makes the next round pay the state-reload cost.
- `iteration_cap` — stop after `iteration_cap` work iterations (default 8). For projects with tight, deterministic rounds.

**If a code review finding is disputed** (AGREE/DISAGREE in `review-protocol.md`): post to `QUESTIONS.md` as an `OPEN` item. Do not silently ignore.

<!-- PROJECT-SPECIFIC: phase-3-extras -->
<!-- Project-specific execution rules. -->
<!-- /PROJECT-SPECIFIC -->

---

## Phase 4 — Close (non-negotiable)

**This phase runs even for aborted rounds.** If you ran out of budget, hit 3 strikes, encountered a stop condition, or are abandoning the priority — you still close the round.

1. **Write the round file.** Create `rounds/round_NNN_<slug>.md` using the template at `.claude/round/round-file-template.md`. `NNN` is zero-padded (001, 002, ...). Slug is a short kebab-case description of the focus.

2. **Update `STATUS.md`:**
   - Append a new row to the Round History table: `NNN | focus | status (done/partial/aborted) | key finding`
   - Update the priority status (TODO → IN PROGRESS → DONE, or add to Dead Ends if abandoned)
   - **REPLACE** the body of `## Current State` with a fresh summary of the new reality (1–3 paragraphs: what works end-to-end, what is in flight, what this round established). Do NOT prepend a new block or accumulate historical state — the round file carries history. STATUS.md must stay bounded; `## Current State` is "what is true right now," not a changelog.
   - Add any new items discovered to `Future Priorities`
   - **Do NOT add accumulating sections** (e.g. "Recent Activity," historical snapshots, or per-round paragraphs outside the defined sections). The canonical STATUS.md structure is: Current State, Active Threads, Future Priorities, Round History, Stop Conditions — nothing else. Accumulated state causes STATUS.md to exceed Phase 1's <30k token budget after ~50 rounds.

3. **Update `KNOWLEDGE_BASE.md`** if this round produced:
   - A new **Validated Decision** (V-entry): phrased as a durable claim with supporting evidence
   - A new **Dead End** (D-entry): an approach confirmed not to work, with reason
   - On audit rounds: write a new **Audit** entry (A-entry) per `review-protocol.md`

4. **Update `QUESTIONS.md`** if this round:
   - Processed any `ANSWERED` items → mark `PROCESSED`
   - Discovered items needing human judgment → add as `OPEN`. Every OPEN item must be phrased as a literal question (title ends in `?`, Short Question block is one plain-English sentence) answerable by a human who has NOT read the round file. Jargon, variable names, file paths, and flag values belong in a trailing `Technical detail` block that future-round workers consume — not the human. See `templates/QUESTIONS.md` for the four-block shape.

5. **Commit.** Single commit per round with targeted `git add` (not `-A`). Before the first `git add`, proactively clear a stale `.git/index.lock` left by an earlier crashed round or aborted hook. This avoids the "Another git process seems to be running" failure that is otherwise recovered reactively from `.claude/round/failure-modes.md` — recognition under context compaction is unreliable, so handle it up front.

   ```bash
   # Stale git index lock cleanup (safe-only: old, no live git, repo not mid-op).
   if [[ -f .git/index.lock ]]; then
     LOCK_AGE=$(( $(date +%s) - $(stat -c %Y .git/index.lock 2>/dev/null || stat -f %m .git/index.lock 2>/dev/null) ))
     if [[ $LOCK_AGE -ge 60 ]] \
        && ! pgrep -u "$USER" -x git >/dev/null 2>&1 \
        && git rev-parse --verify HEAD >/dev/null 2>&1 \
        && [[ ! -d .git/rebase-merge && ! -d .git/rebase-apply && ! -f .git/MERGE_HEAD ]]; then
       echo "Removing stale .git/index.lock (age ${LOCK_AGE}s, no live git process, repo idle)" >&2
       rm -f .git/index.lock
     else
       echo "ERROR: .git/index.lock present and not safely stale — aborting round close." >&2
       echo "  Lock age: ${LOCK_AGE}s (need >=60s to auto-clear)" >&2
       echo "  Investigate: pgrep -u \"\$USER\" -x git; ls -la .git/index.lock; git status" >&2
       exit 5
     fi
   fi

   git add STATUS.md KNOWLEDGE_BASE.md QUESTIONS.md rounds/ <specific changed source files>
   git commit -m "$(cat <<'EOF'
   Round NNN: <focus> — <key finding>

   <one-paragraph summary of what happened, evidence, next steps>

   Co-Authored-By: Claude Opus 4.6 (1M context) <noreply@anthropic.com>
   EOF
   )"
   ```

   **Why the heuristic.** Three signals together distinguish a truly stale lock from a concurrent user operation: (a) age ≥ 60s — active git commands complete in <1s; (b) no live git process for this user — a live holder means it's not stale; (c) repo is not mid-rebase/merge — those can leave `.git/index.lock` legitimately for longer. All three must be true before we remove. On fail, the round aborts cleanly rather than racing; the human can then investigate.

6. **Clear the dashboard active-pick file.** Counterpart to the Phase 2 announcement:
   ```bash
   rm -f /tmp/smlmclustering-active.md
   ```
   The dashboard's status pane shows the active-pick block when the file exists and drops it when the file is gone. Single Bash call; don't use `Edit` or `Write`. Run this before lock release so the pane clears promptly rather than flickering the previous round's pick during the next dispatch's brief pre-lock interval.

7. **Release the lock.** `rm -f "$LOCK"` (the trap also does this on any exit, but release explicitly to signal the dispatcher that the round is cleanly done).

<!-- PROJECT-SPECIFIC: phase-4-extras -->
8. **Periodic @analysis feedback.** After the commit and before releasing the lock, send a terse one-line recap to @analysis and request feedback. This keeps the SMLMAnalysis coordinator informed of interface shifts as SMLMClustering is built out so drift is caught early.

   ```bash
   ~/.claude/skills/agent-comm/scripts/agent-send.sh @analysis \
     "[update] Round NNN (<focus>): <status>. <one-line key finding>. Next: <next-pick>. Feedback welcome on interface/priorities." \
     || true
   ```

   Rules:
   - Keep the message to a single line. Round file, STATUS.md, and KB are the durable record; agent-send is a notification, not a narrative.
   - `|| true` prevents a send failure (network, script error, agent offline) from aborting round close.
   - Do NOT wait for a reply. @analysis's response lands in the inbox and is read in the next round's Phase 1 extras.
   - If the round produced a genuine design question for @analysis, still post it to `QUESTIONS.md` OPEN for Keith — @analysis can weigh in there too, but the human is the decision authority.
<!-- /PROJECT-SPECIFIC -->

---

## Phase 5 — Exit

Print a one-paragraph summary: status (done / partial / aborted / stopped), focus, key finding, pointer to the next priority. This goes to the tmux pane where the dispatcher will read it on window close.

```
--- Round NNN complete ---
Status: <done|partial|aborted|stopped>
Focus: <what this round worked on>
Key finding: <the single most important takeaway>
Next: <what the next round should pick up>
```

Then exit the session (the dispatcher will kill the tmux window when the lock is released).
