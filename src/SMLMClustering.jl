"""
    SMLMClustering

Clustering and spatial-statistic backends for single-molecule localization
microscopy (SMLM) data, operating on `SMLMData.BasicSMLD` emitters.

Three verbs, each dispatched on a concrete config type:

- `cluster` — labeling backends (DBSCAN, HDBSCAN, hierarchical, Voronoi/SR-Tesseler,
  MRF density-regime, point-hysteresis) that write a cluster id onto each
  `emitter.id` (`0` = noise, `1..K` = clusters).
- `cluster_statistics` — read-only spatial statistics (Hopkins clustering tendency,
  Voronoi per-emitter density, local-contrast feature).
- `classify_emitters` — edge / membrane / interior classification.

# Entry point
```julia
(smld_out, info) = cluster(smld, cfg)
```

where `cfg` is a concrete `AbstractClusterConfig` subtype supplied by one of the
backends. See the documentation for the full method catalog.
"""
module SMLMClustering

using SMLMData
using Clustering
using Distances
using DelaunayTriangulation
using NearestNeighbors
using Random

export AbstractClusterConfig, ClusterInfo, cluster,
       DBSCANConfig, HDBSCANConfig, HierarchicalConfig, VoronoiConfig,
       MRFDensityClusterConfig, calibrate_regime_gaussians,
       calibrate_regime_thresholds,
       PointHysteresisConfig,
       AbstractStatisticsConfig, ClusterStatisticsInfo, cluster_statistics,
       HopkinsConfig, VoronoiDensityConfig, LocalContrastFeature

include("types.jl")
include("utils.jl")
include("backends/dbscan.jl")
include("backends/hdbscan.jl")
include("backends/hierarchical.jl")
include("backends/voronoi.jl")
include("backends/hopkins.jl")
include("backends/voronoi_density.jl")
include("backends/local_contrast.jl")
include("backends/mrf_density.jl")
include("backends/point_hysteresis.jl")

include("edge_classify/EdgeClassify.jl")
using .EdgeClassify: classify_emitters,
                     AbstractEdgeClassifyConfig,
                     OuterPolygonConfig, KdeValleyConfig,
                     EdgeClassifyInfo, LoopDiagnostic,
                     in_cell, interior_fraction, method_name, write_edge_artifacts,
                     compute_concavity_metric, ConcavityMetricReport
export classify_emitters,
       AbstractEdgeClassifyConfig,
       OuterPolygonConfig, KdeValleyConfig,
       EdgeClassifyInfo, LoopDiagnostic,
       in_cell, interior_fraction, method_name, write_edge_artifacts,
       compute_concavity_metric, ConcavityMetricReport

end
