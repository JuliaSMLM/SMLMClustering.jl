"""
    EdgeClassify

Edge / membrane / interior classification for 2D SMLM emitter point clouds.

`classify_emitters(in, cfg)` is a top-level verb parallel to the package's
`cluster` / `cluster_statistics`: the concrete config type selects the strategy by
dispatch, and the result is an [`EdgeClassifyInfo`](@ref).

Strategies:
- [`OuterPolygonConfig`](@ref) — multi-K density gate → alpha-shape outer loop →
  point-in-polygon + membrane band.
- [`KdeValleyConfig`](@ref) — validated adaptive dSTORM gate (Gaussian-KDE +
  background/cell valley + footprint fill + enclosure reclass).
"""
module EdgeClassify

using NearestNeighbors
using DelaunayTriangulation
using Statistics
using Dates
import SMLMData

export classify_emitters,
       AbstractEdgeClassifyConfig,
       OuterPolygonConfig, KdeValleyConfig,
       EdgeClassifyInfo, LoopDiagnostic,
       in_cell, interior_fraction, method_name,
       write_edge_artifacts,
       compute_concavity_metric, ConcavityMetricReport

include("configs.jl")
include("info.jl")
include("geometry.jl")
include("gates.jl")
include("diagnostics.jl")
include("classify.jl")
include("io.jl")
include("concavity_metric.jl")

end # module EdgeClassify
