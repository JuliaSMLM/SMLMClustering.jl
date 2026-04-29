# Round Anti-Patterns (Canonical)

Universal guardrails drawn from PSFLearning and SMLMDeepFit evidence. **These apply to every round, in every project, and cannot be overridden.** Project-specific additions go in `.claude/round/project-anti-patterns.md` and are read alongside this file.

If you are about to take an action that looks like an item on this list, stop. Post the situation to `QUESTIONS.md` as `OPEN` and end the round cleanly (Phase 4 still required).

---

## 1. Do not remove or disable tests to make errors go away

If a test fails, the test is giving you information. Removing it discards the information. Skipping it (`@test_skip`, `@test_broken`) is only acceptable when you have already documented WHY the underlying behavior is broken and added an entry to `STATUS.md` Known Issues.

## 2. Do not skip Phase 4 (round close) because the round was short or uneventful

Every round writes a round file, updates `STATUS.md`, and commits. Even a round that accomplished nothing records that it accomplished nothing, and why. The audit trail is the point — skipping close for "easy" rounds destroys the audit trail.

## 3. Do not silently delete KB entries

A Validated Decision or Dead End entry that turned out to be wrong does not get deleted. Mark it `RETIRED` inline with a back-reference to the round file that contradicted it. The history of what was once believed is part of the institutional memory.

## 4. Do not ignore code review findings

When `/review-code` returns findings, each one gets an AGREE or DISAGREE marking (see `review-protocol.md`). Disagreements are posted to `QUESTIONS.md` for the human. Nothing is silently dropped on the floor.

## 5. Do not add dependencies without human approval

Adding to `Project.toml` or `[deps]` requires an `OPEN` question in `QUESTIONS.md`. Dependencies are a long-tail liability; the human decides what the project takes on.

## 6. Do not modify CLAUDE.md, start-round.md, dispatch-round.md, or canonical round files

`CLAUDE.md` encodes project-level rules the human owns. `.claude/commands/start-round.md` and `.claude/commands/dispatch-round.md` are protocol artifacts. `.claude/round/anti-patterns.md`, `failure-modes.md`, `round-file-template.md`, `review-protocol.md` are canonical. If you think one of these files is wrong, post to `QUESTIONS.md`. Editable files are `STATUS.md`, `KNOWLEDGE_BASE.md`, `QUESTIONS.md`, new round files, source code, and `.claude/round/project-*.md`.

## 7. Do not scope-creep

If you discover something worth doing that isn't the current priority — good, write it down. Add it to `Future Priorities` in `STATUS.md` with the full feature trio (`add:`/`test:`/`doc:`). Do NOT pivot the round to chase it. One round, one priority.

## 8. Do not have the orchestrator do tool work directly

(For the orchestrator session — does not apply to the round worker fork.) Every Read / Edit / Bash / Grep call you make in the orchestrator bakes that file content into your context for the rest of the session — even after the file changes on disk. That's pollution. Spawn a fork instead and let it return a short summary. Multi-file reads, dashboard refreshes, queue checks, q-answer processing, "show me priority N" requests — all forks. Trivial single-line edits to a state file where the developer has supplied all needed text inline are OK to do directly.

## 9. Do not commit broken tests as "done"

A priority is not `DONE` until the feature trio's `test:` cases pass. Partial completion with failing tests is `IN PROGRESS` or `BLOCKED`, not `DONE`. The Round History table records partial rounds honestly.

## 10. Do not use external consultation before an in-tree attempt

`/ask-codex` and `/second-opinion` are powerful but expensive. They are allowed only AFTER you have attempted the problem in-tree and failed. Budget: 2 consultations per round. Each recorded in the round file with source and outcome.

## 11. Do not treat warnings as cosmetic

If something is flagged as a warning (lint, type, deprecation), it is not automatically safe to ignore. Warnings are information. Decide: fix it, AGREE, or escalate via QUESTIONS.md — same protocol as review findings.

## 12. Do not improvise the KB audit

The audit procedure is defined concretely in `review-protocol.md`. Follow it step by step. Do not invent a deeper audit because "it feels like more review is needed" — that is ritual inflation.

## 13. Do not post OPEN questions as design dossiers

When adding an item to `QUESTIONS.md` OPEN, the title must end in `?` and the Short Question block must be one plain-English sentence answerable by someone who has not read the round file. Context/proposal/impact dossiers phrased as topics (not questions) put the burden of extracting the actual question onto the human, who often cannot. Jargon — variable names, file paths, flag values, code snippets — belongs in a trailing `Technical detail` block that future-round workers consume, NOT in the question itself. See `templates/QUESTIONS.md` for the four-block shape. If a human has to read your Technical detail to figure out what you're asking, rewrite the Short Question.

## 14. Do not BLOCK, drop scope, skip a priority, or record a Dead End without posting an OPEN question

If a round marks a priority BLOCKED, abandons scope on a priority that was being worked, skips a priority for missing the feature trio, records a new Dead End in `KNOWLEDGE_BASE.md`, or processes an agent message asking for human judgment, an OPEN question in `QUESTIONS.md` is **required** before Phase 4 commit. These are the exact moments where human judgment is load-bearing:

- BLOCKED → "drop scope, wait for external progress, or pivot?"
- Dropped scope → "was the drop right, and what replaces it?"
- Skipped priority (missing trio) → "specify add/test/doc for this priority, or remove it?"
- New Dead End → "pivot to which alternative, or retire the goal?"
- Agent message requesting input → the cross-repo question itself

A BLOCKED priority without a linked OPEN question is a silent dead-letter — the human doesn't see the decision point until they manually audit, which is exactly the audit burden the round system is supposed to remove.

**Exception:** if the identical question is already OPEN from an earlier round (verify by title match against `QUESTIONS.md`), do not post a duplicate. Reference the existing Q in the round file's Next steps section instead.

## 15. Do not add under-specified priorities to STATUS.md

Every new entry in Future Priorities MUST have all three feature-trio legs: `add:`, `test:`, `doc:`. If you can't fully specify the trio (because the work's scope is genuinely unclear), post an OPEN question instead and add the priority only after the question is answered. Adding a vague priority guarantees the next round skips it for missing the trio and posts the same OPEN question — wasted work.
