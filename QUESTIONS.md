# QUESTIONS — SMLMClustering

*Async human-in-the-loop. Rounds post OPEN questions when they need human judgment (disputed review findings, dependency additions, design decisions). Human answers between rounds by editing items into ANSWERED. Next round processes ANSWERED and marks PROCESSED.*

---

## OPEN

<!-- Items awaiting human input. Each item MUST be written for a human who
     has not read any recent round files. Four-block shape, in this order:

### Q1 — <title that ends in a question mark?>

**From round:** NNN

**Short question:** <one plain-English sentence, answerable in isolation — no jargon, no variable names, no file paths>

**What each answer means:**
- **Yes / Option A:** <what changes in the project if the human picks this>
- **No / Option B:** <what changes if they pick this>
- (more options as needed)

**Technical detail:** <free-form context for future-round workers: variable names, file paths, flag values, reference commits, code snippets. This block is NOT read by the human answering the question; it exists so the round that PROCESSES the answer has enough context to act.>

HARD RULE: A human who has not read the round file should be able to answer the Short Question using only the What-each-answer-means block. If they need to read Technical detail to answer, the Short Question is wrong — rewrite it.
-->

(none)

---

## ANSWERED

<!-- Human has responded. Next round will process these and move them to PROCESSED. -->
<!-- Format keeps everything from OPEN and adds:
**Answer:** <human's direction>
**Answered:** YYYY-MM-DD
-->

(none)

---

## PROCESSED

<!-- Archived after a round applied the answer. Keep for audit trail. -->

(none)
