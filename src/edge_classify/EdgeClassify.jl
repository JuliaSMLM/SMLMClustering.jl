"""
    EdgeClassify

Edge / membrane / interior classification for 2D SMLM emitter point clouds.

`classify_emitters(in, cfg)` is a top-level verb parallel to the package's
`cluster` / `cluster_statistics`: the concrete config type selects the strategy by
dispatch, and the result is an [`EdgeClassifyInfo`](@ref).

Strategies:
- [`OuterPolygonConfig`](@ref) — multi-K density gate → alpha-shape outer loop →
  point-in-polygon + membrane band.
- [`GridHybridConfig`](@ref) — outer-polygon + density-grid membrane promotion.
- [`MaskCarveConfig`](@ref) — outer-polygon with a density-mask carve as the
  effective polygon.
- [`KdeValleyConfig`](@ref) — validated adaptive dSTORM gate (Gaussian-KDE +
  background/cell valley + footprint fill + enclosure reclass).

The contract for class labels, fields, and artifacts is documented in
`docs/src/edge_classify_interface_v1.md`.
"""
module EdgeClassify

using NearestNeighbors
using DelaunayTriangulation
using Statistics
using Dates
import SMLMData

export classify_emitters,
       AbstractEdgeClassifyConfig, AbstractPolygonConfig,
       OuterPolygonConfig, GridHybridConfig, MaskCarveConfig, KdeValleyConfig,
       EdgeClassifyInfo, LoopDiagnostic, MaskCarveDiagnostic,
       in_cell, interior_fraction, method_name,
       write_edge_artifacts,
       compute_concavity_metric, ConcavityMetricReport

include("configs.jl")
include("info.jl")
include("geometry.jl")
include("gates.jl")
include("refine.jl")
include("diagnostics.jl")
include("classify.jl")
include("io.jl")
include("concavity_metric.jl")

end # module EdgeClassify
