# Project-Specific Failure Modes — SMLMClustering

Project-local additions to the canonical failure-mode list. Read in addition to `.claude/round/failure-modes.md` when a round encounters unexpected state.

Add project-specific failure modes here as they emerge. Match the canonical format: **Symptom**, **Cause**, **Recovery**.

---

<!-- Example:

## GPU out-of-memory during training step

**Symptom:** `CUDA.CuError: CUDA_ERROR_OUT_OF_MEMORY` during forward pass.

**Cause:** Another process holding GPU memory, or our batch size is too large for the allocated device.

**Recovery:** Check `nvidia-smi`. If another job is running, wait or move to descent. If single-process, reduce batch_size in config and retry.

-->

(add project-specific failure modes here)
