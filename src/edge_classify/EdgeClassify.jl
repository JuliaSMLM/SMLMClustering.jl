"""
    EdgeClassify

Edge / membrane / interior classification for 2D SMLM emitter point clouds.

`classify_emitters(in, cfg)` is a top-level verb parallel to the package's
`cluster` / `cluster_statistics`: the concrete config type selects the strategy by
dispatch, and the result is an [`EdgeClassifyInfo`](@ref).

Strategies:
- [`OuterPolygonConfig`](@ref) — multi-K density gate → alpha-shape outer loop →
  point-in-polygon + membrane band.
- [`KdeValleyConfig`](@ref) — adaptive dSTORM density-valley gate (Gaussian-KDE +
  background/cell valley + footprint fill + enclosure reclass).
"""
module EdgeClassify

using NearestNeighbors
using AdaptivePredicates
using Statistics
using Dates
import SMLMData
import ..SMLMClustering: _point_in_polygon, build_mask, in_region, CellPolygon, MultiCellMask

export classify_emitters,
       AbstractEdgeClassifyConfig,
       OuterPolygonConfig, KdeValleyConfig,
       EdgeClassifyInfo, LoopDiagnostic,
       in_cell, interior_mask, interior_fraction, method_name,
       write_edge_artifacts,
       compute_concavity_metric, ConcavityMetricReport

include("configs.jl")
include("info.jl")
include("delaunay.jl")
include("geometry.jl")
include("gates.jl")
include("diagnostics.jl")
include("classify.jl")
include("io.jl")
include("concavity_metric.jl")

end # module EdgeClassify
