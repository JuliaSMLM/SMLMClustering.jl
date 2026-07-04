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
PrecisionDBSCANConfig
HDBSCANConfig
HierarchicalConfig
VoronoiConfig
MRFDensityClusterConfig
PointHysteresisConfig
calibrate_regime_gaussians
calibrate_regime_thresholds
```

## Precision-DBSCAN primitive

The lower-level, reuse-the-graph primitive behind [`PrecisionDBSCANConfig`](@ref):
build the σ-aware neighbor cache once, then relabel it many times with varying
`σ_eff` / `nsigma` without rebuilding the tree (see the
[Precision DBSCAN](@ref) method page for the walk-through). These names are **public
but not exported** — call them qualified (`SMLMClustering.build_precision_neighbor_graph`,
etc.).

```@docs
PrecisionNeighborGraph
build_precision_neighbor_graph
precision_dbscan_labels
precision_dbscan_labels!
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
in_cell
interior_mask
interior_fraction
```

## Multi-cell mask

```@docs
CellPolygon
MultiCellMask
build_mask
in_region
region_area
```

## Edge-mask report & figures

`compute_edge_report` / `write_edge_report` / `class_codes` are core (no plotting deps).
`plot_edge_report` and `render_classes` are provided by `SMLMClusteringFiguresExt` and
require both `CairoMakie` and `SMLMRender` to be loaded.

```@docs
EdgeReport
compute_edge_report
write_edge_report
class_codes
```

## Advanced & diagnostics

```@docs
method_name
write_edge_artifacts
compute_concavity_metric
ConcavityMetricReport
LoopDiagnostic
```
