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

## 6. Do not modify CLAUDE.md, start-round.md, or canonical round files

`CLAUDE.md` encodes project-level rules the human owns. `.claude/commands/start-round.md` is the protocol artifact. `.claude/round/anti-patterns.md`, `failure-modes.md`, `round-file-template.md`, `review-protocol.md` are canonical. If you think one of these files is wrong, post to `QUESTIONS.md`. Editable files are `STATUS.md`, `KNOWLEDGE_BASE.md`, `QUESTIONS.md`, new round files, source code, and `.claude/round/project-*.md`.

## 7. Do not scope-creep

If you discover something worth doing that isn't the current priority — good, write it down. Add it to `Future Priorities` in `STATUS.md`. Do NOT pivot the round to chase it. One round, one priority.

## 8. Do not bypass the lock file

The lock is how `/loop /dispatch-round` prevents overlapping rounds. Manually deleting the lock file while a round is running will corrupt the system. The only legitimate reason to delete the lock is a confirmed crash (tmux window dead but lock still there — see `failure-modes.md`).

## 9. Do not commit broken tests as "done"

A priority is not `DONE` until tests covering it pass. Partial completion with failing tests is `IN PROGRESS` or `BLOCKED`, not `DONE`. The Round History table records partial rounds honestly.

## 10. Do not use external consultation before an in-tree attempt

`/ask-codex` and `/second-opinion` are powerful but expensive. They are allowed only AFTER you have attempted the problem in-tree and failed. Budget: 2 consultations per round. Each recorded in the round file with source and outcome.

## 11. Do not treat warnings as cosmetic

If something is flagged as a warning (lint, type, deprecation), it is not automatically safe to ignore. Warnings are information. Decide: fix it, AGREE, or escalate via QUESTIONS.md — same protocol as review findings.

## 12. Do not improvise the KB audit

The audit procedure is defined concretely in `review-protocol.md`. Follow it step by step. Do not invent a deeper audit because "it feels like more review is needed" — that is ritual inflation.

## 13. Do not post OPEN questions as design dossiers

When adding an item to `QUESTIONS.md` OPEN, the title must end in `?` and the Short Question block must be one plain-English sentence answerable by someone who has not read the round file. Context/proposal/impact dossiers phrased as topics (not questions) put the burden of extracting the actual question onto the human, who often cannot. Jargon — variable names, file paths, flag values, code snippets — belongs in a trailing `Technical detail` block that future-round workers consume, NOT in the question itself. See `templates/QUESTIONS.md` for the four-block shape. If a human has to read your Technical detail to figure out what you're asking, rewrite the Short Question.
