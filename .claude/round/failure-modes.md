# Round Failure Modes (Canonical)

Known failure modes and their recovery procedures. Consult this file when a round encounters unexpected state. Project-specific additions go in `.claude/round/project-failure-modes.md`.

---

## Orphaned lock (lock file exists, no round running)

**Symptom:** `/dispatch-round` reports lock held, but no `<project>-round` tmux window exists and no round process is alive.

**Cause:** Previous round crashed (process killed, tmux window killed externally, system reboot) before Phase 4 released the lock.

**Recovery:**
- If `lock_strategy: pid_based` — the dispatcher should reclaim automatically (check PID in lock file against `kill -0`).
- If `lock_strategy: manual_reclamation` — inspect the lock file (`cat /tmp/<project>-round.lock`), confirm no round is active, then delete manually: `rm /tmp/<project>-round.lock`.
- If `STATUS.md` has an incomplete Round History entry for the crashed round, add a note that the round aborted without closing cleanly. No retroactive round file.

---

## Orphaned tmux window (window exists, no lock)

**Symptom:** `<project>-round` tmux window exists but the lock file is gone, or the window is idle at a shell prompt.

**Cause:** Round closed cleanly (Phase 4 released lock) but the dispatcher monitor loop didn't run the cleanup (dispatcher killed, network drop, etc.).

**Recovery:** `tmux kill-window -t <project>-round`. The next `/dispatch-round` invocation will create a fresh window.

---

## Git index lock (`.git/index.lock`)

**Symptom:** Phase 4 commit fails with "Another git process seems to be running."

**Cause:** A previous round crashed mid-commit, leaving `.git/index.lock`. Or a concurrent git operation outside the round.

**Recovery:** Confirm no git process is running (`ps aux | grep git`), then `rm .git/index.lock` and retry the commit. If git operations are happening outside the round, the round's commit should wait rather than racing.

---

## Round ran without committing

**Symptom:** Round file exists in `rounds/` but `git log` shows no corresponding `Round NNN:` commit. `STATUS.md` may or may not be updated.

**Cause:** Phase 4 hit an error between writing files and committing. Very rare.

**Recovery:** Inspect `git status`. The round file and `STATUS.md` changes should be visible as uncommitted changes. Make the commit manually with the same message format. If `STATUS.md` wasn't updated, add the Round History row manually before committing.

---

## Tests pass individually but fail in the full suite

**Symptom:** A priority's targeted test passes in isolation, but `Pkg.test()` fails.

**Cause:** Test pollution — global state, load order, file system assumptions.

**Recovery:** This is a legitimate round-3 strike material. Do not mark the priority `DONE` until the full suite passes. If 3 attempts fail, document as a Dead End: the specific interaction that causes the pollution, a workaround if any, and escalate via `QUESTIONS.md`.

---

## Opaque errors from Reactant / Enzyme / compiler-heavy libraries

**Symptom:** Stack trace lands inside a compilation step, no clear actionable line in user code.

**Cause:** Compiler-internal failures often have fixes outside the code you're working on (version mismatch, cached IR, upstream bug).

**Recovery:** First check versions and try clearing compile caches. If still opaque after one in-tree attempt, this is explicit external-consultation material: use `/ask-codex` or `/second-opinion` (within the 2/round budget). Record the consulted source and outcome in the round file.

---

## 3-strike Dead End reached

**Symptom:** You have attempted an obstacle 3 times with different approaches, all failed.

**Cause:** Working as designed. This is the protocol's escape valve against infinite rabbit holes.

**Recovery:** Stop attempting. Write a new Dead End (D-entry) in `KNOWLEDGE_BASE.md` describing the obstacle, the three approaches tried, and why each failed. Update the current priority in `STATUS.md` to `BLOCKED` with a reference to the new D-entry. Post an `OPEN` question to `QUESTIONS.md` asking the human for direction. Close the round normally.

---

## Stop Condition met in Phase 2

**Symptom:** Phase 2 detects a Stop Condition in `STATUS.md`.

**Cause:** Human has set a stop condition intentionally (pause autonomous work). Or a protocol condition triggered (e.g., all priorities DONE, test suite failing project-wide).

**Recovery:** Halt with `ROUND STOPPED: <which condition>`. Still run Phase 4 with an aborted round file recording the stop. The lock releases normally. The next `/dispatch-round` will check conditions again and skip again until the human removes the condition.

---

## Dispatcher orphan detection fires

**Symptom:** Dispatch script reports "orphaned lock, window missing — manual intervention required" and refuses to proceed.

**Cause:** Lock held, tmux window gone. Round crashed.

**Recovery:** Follow the "orphaned lock" procedure above. Delete the lock file after confirming no round is running.
