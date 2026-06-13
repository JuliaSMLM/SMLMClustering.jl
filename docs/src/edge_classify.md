# Edge / Membrane Classification

Classify each 2D SMLM emitter as `:outside`, `:membrane`, or `:interior` — the
off-cell background, the cell-boundary band, and the cell interior. The verb
`classify_emitters` is a peer of the package's `cluster` / `cluster_statistics`:
the **concrete config type selects the strategy by dispatch**, and the result is an
`EdgeClassifyInfo`.

## API

```julia
smld, info = classify_emitters(smld::BasicSMLD, cfg::AbstractEdgeClassifyConfig)
info        = classify_emitters(x_um, y_um, cfg::AbstractEdgeClassifyConfig; fov_um)
```

- **SMLD form** (pipeline-facing, `(out, Info)` convention): returns the smld with
  the primary class mirrored into `smld.metadata["edge_classify_class"]`
  (`Vector{String}`), plus the `EdgeClassifyInfo`. `fov_um` is taken from
  `smld.camera.pixel_edges_x/y`.
- **Coordinate form** (computational core): returns the `info` directly.
  `fov_um = (xmin, xmax, ymin, ymax)` in µm (accepts any `Real` tuple).

Artifact writing is a separate step (compute and IO are decoupled):

```julia
write_edge_artifacts(leaf, info, x_um, y_um; condition, cell)
```

## Strategies (configs)

Each config is a `<: AbstractEdgeClassifyConfig` (sibling of `AbstractClusterConfig`)
holding only its own parameters; fields are lowercase, validated at dispatch entry.

### `OuterPolygonConfig`

FOV-reflection → multi-K k-NN density gate → alpha-shape outer loop →
point-in-polygon + a `membrane_nm` band.

| field | default | unit | meaning |
|---|---|---|---|
| `alpha_nm` | 300 | nm | alpha-shape circumradius |
| `membrane_nm` | 100 | nm | membrane band width inboard of the outer polygon |
| `reflect_radius_nm` | 1500 | nm | mirror band inboard of truncated FOV sides |
| `fov_trunc_tol_nm` | 150 | nm | FOV-truncation detection tolerance |
| `k_list` | `(16, 128)` | — | k-NN K values for the multi-K density gate (intersection) |
| `rho_k_thresh` | 200 | µm⁻² | per-K density gate threshold |

### `KdeValleyConfig`

Validated adaptive dSTORM gate. Gaussian-KDE density on the **original** cloud →
background/cell valley threshold → footprint fill → the outer-polygon geometry on
the footprint subset → ray-cast **enclosure** reclass folding enclosed background
into `:interior`. Per-FOV adaptive — no per-cell density tuning. `alpha_nm = 600`
is the validated value (the type carries it; no factory needed).

| field | default | unit | meaning |
|---|---|---|---|
| `alpha_nm` | 600 | nm | alpha-shape circumradius (validated) |
| `membrane_nm` | 100 | nm | membrane band width |
| `reflect_radius_nm` | 1500 | nm | mirror band inboard of truncated sides |
| `fov_trunc_tol_nm` | 150 | nm | FOV-truncation tolerance |
| `sigma_nm` | 150 | nm | Gaussian-KDE bandwidth σ |
| `rmax_sigma` | 3.0 | — | KDE range-query cutoff in units of σ |
| `valley_nbins` | 140 | — | log-density histogram bins for the valley threshold |
| `valley_floorfrac` | 0.05 | — | left-base cutoff as a fraction of the cell-mode peak |
| `valley_smooth` | 4 | bins | ± window for histogram smoothing |
| `footprint_bin_um` | 0.2 | µm | raster bin for the footprint fill |
| `footprint_closing_px` | 3 | px | morphological closing radius (seal thin necks) |
| `enclosure_bin_um` | 0.2 | µm | raster bin for the 8-ray enclosure reclass |
| `enclosure_min_hits` | 6 | of 8 | min rays hitting cell tissue to fold a point into `:interior` |

## Result — `EdgeClassifyInfo`

`EdgeClassifyInfo{C} <: SMLMData.AbstractSMLMInfo`. Key fields:

| field | type | meaning |
|---|---|---|
| `class` | `Vector{Symbol}` | authoritative per-emitter answer: `:outside` / `:membrane` / `:interior` |
| `inside_outer` | `BitVector` | **geometric** containment in the alpha outer loop |
| `dist_to_outer_um` | `Vector{Float64}` | distance to the outer polygon; `NaN` when not inside |
| `outer_polygon`, `loops` | polygons | the alpha outer loop + all loops |
| `loop_diagnostics` | `Vector{LoopDiagnostic}` | per-loop diagnostics |
| `config` | `C` | the concrete config that ran (provenance) |
| `fov_um`, `truncated_sides`, `n_reflected`, `runtime_s` | — | run metadata |
| `n_outside`, `n_membrane`, `n_interior` | `Int` | class counts |

Accessors: `in_cell(info)` = `info.class .!= :outside` (topological membership);
`interior_fraction(info)`.

**Class semantics.** `class` is the canonical answer — **filter on `class`, never on
`inside_outer`**. For `OuterPolygonConfig`, `inside_outer` is geometric for every
emitter and `in_cell(info) == inside_outer`. For `KdeValleyConfig`, the geometry is
computed on the footprint subset (off-footprint emitters carry `inside_outer =
false` / `dist = NaN`) and the enclosure stage folds enclosed background into
`class == :interior` while leaving `inside_outer` geometric — so:

- `in_cell ⊇ inside_outer`,
- the **enclosure-recovered set** is exactly `class == :interior && inside_outer == false`
  (those carry `dist_to_outer_um == NaN`),
- `membrane` is always the band around the geometric outer polygon.

## Class invariants

- `:outside ∪ :membrane ∪ :interior` partitions the input set (no nulls, no overlaps).
- Order matches the input.
- The SMLD form mirrors `class` (as `String`s) into `smld.metadata["edge_classify_class"]`.

## Artifacts (`write_edge_artifacts`)

Written under `<out_dir>/<condition>/<cell>/`. Schemas are stamped in headers /
`manifest.json`; per-config params are serialized via the `to_dict` trait (only the
fields that actually ran are recorded).

| file | schema | contents |
|---|---|---|
| `classified.tsv` | 2 | `emitter_id, x_um, y_um, class, inside_outer, in_cell, dist_to_outer_um` |
| `polygon_loops.tsv` | 1 | all alpha-shape loops (`loop_id, vertex_id, x_um, y_um`) |
| `loop_diagnostics.csv` | 2 | per-loop diagnostics |
| `params.json` | 2 | git provenance, fov, truncation, `params` (method-specific via `to_dict`), runtime |
| `manifest.json` | 1 | artifact index + schema versions |

`params.json` carries `params["METHOD"]` = `method_name(cfg)` (`"outer_polygon"` /
`"kde_valley"`) as a write-only provenance label.

## Concavity metric (diagnostic)

```julia
report = compute_concavity_metric(info, x_um, y_um; buffer_um=2.0, ...)
```

Flags `:interior` emitters that sit in deep concave bays the alpha-shape bridged
across (boundary-proximal, high directional asymmetry, low local density),
stratified by whether the nearest outer segment is inside the FOV or straddles its
edge. Diagnostic only — does not change `class`. Returns a `ConcavityMetricReport`.
