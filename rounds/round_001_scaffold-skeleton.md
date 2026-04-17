# Round 001 — Scaffold package skeleton

**Date:** 2026-04-17
**Status:** done
**Priority worked on:** Priority 1 — Scaffold package skeleton (SMLMData dep + `AbstractClusterConfig` + `ClusterInfo` + `cluster` entry point)

## Hypothesis

The design already agreed with @analysis and captured in `KNOWLEDGE_BASE.md` (V1–V4) is concrete enough to implement as a standalone interface package — abstract config type, result type, and entry-point generic — with no backend wired yet. If that is true, the four backend rounds that follow should only need to add a concrete `*Config <: AbstractClusterConfig` and a single `cluster(smld, cfg::...)` method each, without revisiting the interface.

## What was attempted

- **Project.toml** — added `SMLMData = "5488f106-..."` to `[deps]` and `SMLMData = "0.7"` to `[compat]`. Matches the pattern in `SMLMFrameConnection/Project.toml`.
- **src/types.jl** (new) — defined `abstract type AbstractClusterConfig <: SMLMData.AbstractSMLMConfig`, `struct ClusterInfo <: SMLMData.AbstractSMLMInfo` with the seven fields agreed in STATUS Priority 1 (`n_locs_in, n_clustered, n_noise, n_clusters, cluster_sizes::Vector{Int}, algorithm::Symbol, elapsed_s::Float64`), and a fallback `cluster(smld::BasicSMLD, cfg::AbstractClusterConfig)` method that errors with a pointer to the concrete backends.
- **src/SMLMClustering.jl** — replaced the generator stub with `using SMLMData`, `include("types.jl")`, and `export AbstractClusterConfig, ClusterInfo, cluster`.
- **test/runtests.jl** — three testsets: subtyping (`AbstractClusterConfig <: AbstractSMLMConfig`, `ClusterInfo <: AbstractSMLMInfo`), `ClusterInfo` field construction with a consistency check (`sum(cluster_sizes) == n_clustered`), and an `@test_throws ErrorException` against the abstract `cluster` fallback using a top-level dummy subtype.
- **Dep resolution** — `Pkg.develop(path=...)` pointed at the local `SMLMData` checkout, then `Pkg.resolve/instantiate/test`.

## What worked

- Full test suite passes: 11/11 in 0.1 s (`Pkg.test`).
- Package precompiles cleanly against the in-tree `SMLMData v0.7.0`.
- Test suite exercises the exact fallback path future backends will override, so regressions there (e.g. accidentally making the abstract fallback permissive) will be caught.

## What failed

Nothing material. One minor correction mid-round: first draft of `runtests.jl` defined the dummy config subtype inside a `@testset` block, which Julia ≥1.9 rejects (top-level-only). Moved to file scope before the first testset.

## Files changed

```
Project.toml                                        +3 -0
Manifest.toml                                       regenerated
src/SMLMClustering.jl                               +20 -3
src/types.jl                                        +76 (new)
test/runtests.jl                                    +34 -3
rounds/round_001_scaffold-skeleton.md               new
STATUS.md                                           Current State + priority statuses + Round History
```

## Confidence

High — full test suite passes, and every interface claim in STATUS/KB V1–V4 (`AbstractClusterConfig` subtyping, `ClusterInfo` shape, tuple-returning entry point) has a direct assertion in `test/runtests.jl`. The only pieces not exercised are backend-specific, which is correct for this round.

## External consultations (if any)

None.

## Next steps

Priority 2: implement `DBSCANConfig <: AbstractClusterConfig` + `cluster(smld, cfg::DBSCANConfig)`. Add `Clustering.jl` to `[deps]`. Key design points already decided: per-dataset loop gated by `cfg.per_dataset`, label emitters in-place (remember `Emitter2DFit` / `Emitter3DFit` are `mutable`), honor `cfg.remove_unclustered`, populate `ClusterInfo` from the resulting labels. Write a synthetic-blob test with known ground-truth cluster count.

## Questions posted

None.
