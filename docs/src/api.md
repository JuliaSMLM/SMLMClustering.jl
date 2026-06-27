```@meta
CurrentModule = SMLMClustering
```

# API Reference

Complete reference for the exported API. The user-facing verbs and configuration
types lead; diagnostics and lower-level helpers follow under
[Advanced & diagnostics](@ref).

```@index
```

## Core — labeling

```@docs
cluster
AbstractClusterConfig
ClusterInfo
DBSCANConfig
HDBSCANConfig
HierarchicalConfig
VoronoiConfig
MRFDensityClusterConfig
PointHysteresisConfig
calibrate_regime_gaussians
calibrate_regime_thresholds
```

## Spatial statistics

```@docs
cluster_statistics
AbstractStatisticsConfig
ClusterStatisticsInfo
HopkinsConfig
VoronoiDensityConfig
LocalContrastFeature
```

## Edge classification

```@docs
classify_emitters
AbstractEdgeClassifyConfig
OuterPolygonConfig
KdeValleyConfig
EdgeClassifyInfo
```

## Advanced & diagnostics

```@docs
in_cell
interior_fraction
method_name
write_edge_artifacts
compute_concavity_metric
ConcavityMetricReport
LoopDiagnostic
```
