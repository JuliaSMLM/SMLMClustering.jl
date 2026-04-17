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

export AbstractClusterConfig, ClusterInfo, cluster

include("types.jl")

end
