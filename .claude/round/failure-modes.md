# Round Failure Modes (Canonical)

Known failure modes and their recovery procedures. Consult this file when a round encounters unexpected state. Project-specific additions go in `.claude/round/project-failure-modes.md`.

---

## Round fork interrupted mid-execution

**Symptom:** The round fork was killed (orchestrator interrupted, Agent error, network issue) before reaching Phase 4. Files may be partially edited; no commit landed; `/tmp/<slug>-active.md` may still exist with the previous round's pick.

**Cause:** Orchestrator interrupt, transient API error, fork timeout, machine sleep.

**Recovery:**
- The orchestrator's next dispatch sees no in-flight Agent call → spawns a fresh round fork.
- The new fork's Phase 0 truncates `/tmp/<slug>-worker-progress.md` (overwriting the interrupted run's log).
- The new fork's Phase 2 picks the priority again. If the previous fork made source-file edits but didn't commit, those changes appear as uncommitted in `git status` — the new fork should `git restore` them before starting (they're partial work the previous round abandoned), unless the developer indicates otherwise.
- Stale `/tmp/<slug>-active.md` from the killed round is overwritten in the new round's Phase 2.

If the killed round committed but didn't update STATUS.md (very rare — Phase 4 ordering is STATUS.md → commit), inspect `git log -1` and reconcile manually.

---

## Pre-flight fork returns inconsistent state

**Symptom:** Orchestrator's pre-flight returns garbled output, missing fields, or incoherent counts. The orchestrator can't tell whether the queue is empty.

**Cause:** Transient model error in the small fork, or the project's STATUS.md / QUESTIONS.md has a structural problem (missing section header, malformed Future Priorities).

**Recovery:**
- Retry the pre-flight fork once.
- If still garbled, post to the developer: "Pre-flight returned inconsistent state — likely a STATUS.md structural issue. Skipping this dispatch."
- Skip the round; do NOT terminate the loop. Next /loop tick retries from scratch.
- If two consecutive pre-flights fail, terminate the loop and ask the developer to inspect STATUS.md and QUESTIONS.md.

---

## Orchestrator compacted mid-round

**Symptom:** Orchestrator's context auto-compacted while the round fork was running. The orchestrator's chat summary may have lost discussion details.

**Cause:** Long-lived orchestrator session crossed the auto-compact threshold.

**Recovery:** No action needed. The fork's return arrives at the (now-compacted) orchestrator and the digest synthesis still works — the key facts are in the fork's return value, not in the orchestrator's lost discussion. The orchestrator's discussion narrative is degraded for future rounds but not destroyed; rebuilding via on-disk artifacts (STATUS.md, recent round files, KNOWLEDGE_BASE.md) is automatic on each round's Phase 1 read.

To pre-empt: when the orchestrator session has been running ≥30 rounds, periodically `/clear` after archiving any unique discussion to a project notes file. This is a manual hygiene step, not automated.

---

## Git index lock (`.git/index.lock`)

**Symptom:** Phase 4 commit fails with "Another git process seems to be running."

**Cause:** A previous round crashed mid-commit, leaving `.git/index.lock`. Or a concurrent git operation outside the round.

**Recovery:** Phase 4 has a built-in safe-stale check (age ≥ 60s, no live git process for this user, repo not mid-rebase/merge). If those conditions hold, it auto-clears. Otherwise it aborts the round close with a diagnostic. Investigate manually: `pgrep -u "$USER" -x git`, `ls -la .git/index.lock`, `git status`. If genuinely stale, `rm .git/index.lock` and re-run.

---

## Round ran without committing

**Symptom:** Round file exists in `rounds/` but `git log` shows no corresponding `Round NNN:` commit. `STATUS.md` may or may not be updated.

**Cause:** Phase 4 hit an error between writing files and committing. Very rare.

**Recovery:** Inspect `git status`. The round file and `STATUS.md` changes should be visible as uncommitted changes. Make the commit manually with the same message format. If `STATUS.md` wasn't updated, add the Round History row manually before committing.

---

## Tests pass individually but fail in the full suite

**Symptom:** A priority's targeted test passes in isolation, but `Pkg.test()` fails.

**Cause:** Test pollution — global state, load order, file system assumptions.

**Recovery:** This is legitimate 3-strike material. Do not mark the priority `DONE` until the full suite passes. If 3 attempts fail, document as a Dead End: the specific interaction that causes the pollution, a workaround if any, and escalate via `QUESTIONS.md`.

---

## Opaque errors from Reactant / Enzyme / compiler-heavy libraries

**Symptom:** Stack trace lands inside a compilation step, no clear actionable line in user code.

**Cause:** Compiler-internal failures often have fixes outside the code you're working on (version mismatch, cached IR, upstream bug).

**Recovery:** First check versions and try clearing compile caches. If still opaque after one in-tree attempt, this is explicit external-consultation material: use `/ask-codex` or `/second-opinion` (within the 2/round budget). Record the consulted source and outcome in the round file.

---

## 3-strike Dead End reached

**Symptom:** You have attempted an obstacle 3 times with different approaches, all failed.

**Cause:** Working as designed. This is the protocol's escape valve against infinite rabbit holes.

**Recovery:** Stop attempting. Write a new Dead End (D-entry) in `KNOWLEDGE_BASE.md` describing the obstacle, the three approaches tried, and why each failed. Update the current priority in `STATUS.md` to `BLOCKED` with a reference to the new D-entry. Post an `OPEN` question to `QUESTIONS.md` asking the human for direction. **If you ran as Sonnet**, additionally tag the priority `[needs_opus]` in STATUS.md so the next round's pre-flight escalates. Close the round normally.

---

## Stop Condition met in Phase 2

**Symptom:** Phase 2 detects a Stop Condition in `STATUS.md`.

**Cause:** Human has set a stop condition intentionally (pause autonomous work). Or a protocol condition triggered (e.g., test suite failing project-wide).

**Recovery:** Halt with status `stopped`. Still run Phase 4 with an aborted round file recording the stop. The next `/dispatch-round` pre-flight detects the same stop condition and skips the round (without terminating /loop — the orchestrator keeps polling for the condition to clear). The user removes the condition by editing STATUS.md.

---

## Pre-flight returns ACTIONABLE=0 but UNDER_SPECIFIED>0

**Symptom:** STATUS.md has Future Priorities entries, but none have the `add:`/`test:`/`doc:` trio. The pre-flight reports 0 actionable, N under-specified.

**Cause:** Priorities seeded by a human who didn't know about the trio shape, or migrated from an older STATUS.md format.

**Recovery:** Orchestrator (per `dispatch-round.md` Step 2) spawns a small fork that posts an OPEN question per under-specified priority asking the developer to specify the trio. /loop continues — once the developer answers, the next round processes the answer and re-shapes the priority. The orchestrator does NOT terminate /loop here; the queue still has work pending human input.

---

## Auto-stop fired but the developer has work pending

**Symptom:** /loop terminated with "Queue empty" but the developer thought there were priorities to work on.

**Cause:** The priorities in Future Priorities were under-specified (missing trio) AND the OPEN question count happened to be zero at that moment (e.g. all OPEN became ANSWERED in the previous round but Phase 4 didn't add new OPEN items because the priorities were under-specified for a different reason).

**Recovery:** Edit STATUS.md to add a properly-shaped feature-trio priority. Re-arm with `/loop 15m /dispatch-round`. The next pre-flight will see ACTIONABLE>0 and proceed.

To pre-empt: when seeding STATUS.md initially, always use the trio shape — see `templates/STATUS.md` for the canonical format.
