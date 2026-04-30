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
       AbstractStatisticsConfig, ClusterStatisticsInfo, cluster_statistics,
       HopkinsConfig, VoronoiDensityConfig

include("types.jl")
include("utils.jl")
include("backends/dbscan.jl")
include("backends/hdbscan.jl")
include("backends/hierarchical.jl")
include("backends/voronoi.jl")
include("backends/hopkins.jl")
include("backends/voronoi_density.jl")
include("backends/mrf_density.jl")

end
