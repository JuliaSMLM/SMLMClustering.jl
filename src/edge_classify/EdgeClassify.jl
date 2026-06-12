"""
    EdgeClassify

Edge / membrane / interior classification for 2D SMLM emitter point
clouds. The v1 classifier is the **outer-polygon** classifier:

1. Detect FOV-truncated camera-frame sides.
2. Mirror emitters across truncated sides (FOV-augmented set).
3. Multi-K density gate on augmented set (`K_LIST`, `RHO_K_THRESH`)
   yielding a tissue mask.
4. Alpha-shape (`ALPHA_NM`) on tissue points → boundary loops, sorted by
   `abs(area)` descending.
5. For each ORIGINAL emitter, classify by point-in-polygon vs the outer
   loop (`loop_id == 1`) plus a `MEMBRANE_NM` band: `outside`,
   `membrane`, or `interior`.

Interior loops (`loop_id >= 2`) are diagnostic only in v1.

# Public API

- [`classify_emitters`](@ref) — coordinate-based core function and SMLD
  adapter.
- [`EdgeClassifyConfig`](@ref) — parameter struct (defaults provisional).
- [`EdgeClassificationResult`](@ref) — result type.
- [`LoopDiagnostic`](@ref) — per-loop diagnostic record.

The contract for inputs, outputs, filenames, schemas, and class
invariants is documented in `docs/src/edge_classify_interface_v1.md`.
"""
module EdgeClassify

using NearestNeighbors
using DelaunayTriangulation
using Statistics
using Dates

export classify_emitters,
       EdgeClassifyConfig, EdgeClassifyParams,
       EdgeClassificationResult, LoopDiagnostic,
       MaskCarveDiagnostic, kde_valley_params,
       compute_concavity_metric, ConcavityMetricReport

include("types.jl")
include("geometry.jl")
include("grid_hybrid.jl")
include("mask_carve.jl")
include("kde_valley.jl")
include("diagnostics.jl")
include("io.jl")
include("classify.jl")
include("smld_adapter.jl")
include("concavity_metric.jl")

end # module EdgeClassify
