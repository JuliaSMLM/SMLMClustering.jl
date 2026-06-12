"""
    SMLMClustering

Clustering backends for single-molecule localization microscopy data.

Provides a common interface over four clustering algorithms — DBSCAN,
HDBSCAN, Voronoi tessellation, and hierarchical — operating on
`SMLMData.BasicSMLD` emitters. Each backend writes cluster labels back
to `emitter.id` (`0` = noise, `1..K` = cluster).

# Entry point
```julia
(smld_out, info) = cluster(smld, cfg)
```

where `cfg` is a concrete `AbstractClusterConfig` subtype supplied by
one of the backends.
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
                     AbstractEdgeClassifyConfig, AbstractPolygonConfig,
                     OuterPolygonConfig, GridHybridConfig, MaskCarveConfig, KdeValleyConfig,
                     EdgeClassifyInfo, LoopDiagnostic, MaskCarveDiagnostic,
                     in_cell, interior_fraction, method_name, write_edge_artifacts,
                     compute_concavity_metric, ConcavityMetricReport
export classify_emitters,
       AbstractEdgeClassifyConfig, AbstractPolygonConfig,
       OuterPolygonConfig, GridHybridConfig, MaskCarveConfig, KdeValleyConfig,
       EdgeClassifyInfo, LoopDiagnostic, MaskCarveDiagnostic,
       in_cell, interior_fraction, method_name, write_edge_artifacts,
       compute_concavity_metric, ConcavityMetricReport

end
