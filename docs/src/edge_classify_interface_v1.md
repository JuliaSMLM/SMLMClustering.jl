# Edge / Membrane Classification — Interface v1 (locked)

Status: **LOCKED — signed off by @codex-cluster**. This document defines
the v1 contract for the per-emitter edge/membrane/interior classification
pipeline so @genmab and @analysis can build around it. Implementation may
proceed.

Revision history:
- v0.1 — initial draft
- v0.2 — aligned to @codex-cluster directives (coordinate-based core API,
  nested output layout, manifest.json, outer-only class decision,
  dist_to_outer_um column, loop_diagnostics includes outer loop, --out
  authoritative)
- **v1 (this revision)** — final tighten-up: `type` → `heuristic_type`
  (placed last in `loop_diagnostics.csv`); `fov_um` order pinned to
  `(xmin_um, xmax_um, ymin_um, ymax_um)` with validation; params.toml
  uses uppercase keys, unknown keys error

---

## 1. Module & entry points

Module path: `SMLMClustering.EdgeClassify`.

### 1a. Stable core API (coordinate / FOV based)

```julia
using SMLMClustering.EdgeClassify

result = classify_emitters(
    x_um::Vector{Float64},
    y_um::Vector{Float64};
    fov_um::NTuple{4,Float64},          # (xmin, xmax, ymin, ymax)
    params::EdgeClassifyParams = EdgeClassifyParams(),
    out_dir::Union{Nothing,String} = nothing,
    condition::Union{Nothing,String} = nothing,
    cell::Union{Nothing,String} = nothing,
    write_artifacts::Bool = false,
    write_renders::Bool = false,
)::EdgeClassificationResult
```

If `write_artifacts=true`, `out_dir`, `condition`, and `cell` are all
required (the function errors otherwise) — they determine the output path
in §3.

### 1b. SMLD adapter

```julia
result = classify_emitters(
    smld;                               # SMLMData.SMLD or path to JLD2
    params = EdgeClassifyParams(),
    out_dir = nothing,
    condition = nothing,
    cell = nothing,
    write_artifacts = false,
    write_renders = false,
)::EdgeClassificationResult
```

The adapter extracts `x_um`, `y_um` from `smld.emitters[].x/.y`, and
`fov_um` from `smld.camera.pixel_edges_x[1/end]`, `pixel_edges_y[1/end]`,
then calls the core API. No other behavior differences.

### 1c. CLI script

```
julia --project=dev/scripts dev/scripts/edge_classify.jl \
    --smld      <path-to-smld_bagol.jld2>                \
    --condition <COND>                                    \
    --cell      <CELL>                                    \
    --out       <out_dir>                                 \
    [--params   <params.toml>]                            \
    [--renders]
```

The CLI **always** writes artifacts. `--out` is **authoritative**: outputs
go strictly under that path (no fallback to `dev/scripts/output` or any
script-relative working directory).

`params.toml` accepts these uppercase keys (any subset; missing keys use
defaults; unknown keys error): `K_LIST`, `RHO_K_THRESH`, `ALPHA_NM`,
`REFLECT_RADIUS_NM`, `MEMBRANE_NM`, `FOV_TRUNC_TOL_NM`, `METHOD`,
`GRID_PX_NM`, `GRID_SMOOTH_NM`, `GRID_MASK_Q`,
`GRID_MASK_PEAK_FRAC`, `GRID_OUTER_BUFFER_NM`,
`CONCAVITY_METRIC_BUFFER_NM`. The fully resolved parameter set is
written to `params.json`.

### 1d. Concavity metric (Stage 1, this branch)

```julia
report = compute_concavity_metric(result, x_um, y_um;
    buffer_um = result.params_used.CONCAVITY_METRIC_BUFFER_NM/1000,
    asym_R_nm = 1000.0, asym_gate = 0.20, rho_lo = 200.0,
    intracellular_dense_threshold = 0.8)::ConcavityMetricReport
```

Boundary-proximal, loop-aware concavity-error metric for the v1
outer-polygon classifier. Identifies v1 `interior` emitters that are
within `buffer_um` of the outer polygon, have high directional asymmetry
at long radius, and low local density — the chord-vertex signature of a
deep concave bay the alpha-shape bridged across.

Excludes emitters inside intracellular-void loops (raw-metric: non-outer
loops with `frac_in_fov == 1.0` AND `frac_dense >=
intracellular_dense_threshold`), so nuclei / sparse intracellular regions
are not counted as membrane concavity errors. Stratifies suspects by
whether the nearest outer-polygon segment is interior-FOV (Approach B
target — vertex snap to asym ridge) or FOV-edge (Approach C target —
FOV-crossing loops as auxiliary boundary).

The metric does **not** modify classification; it is read-only post-hoc.

---

## 2. Inputs

Callers pass **originals only**. The stable v1 contract does NOT accept
pre-reflected emitters.

| Field | Source (SMLD adapter) | Units | Notes |
|-------|-----------------------|-------|-------|
| `x_um`, `y_um` | `smld.emitters[i].x/.y` | µm | per-emitter coordinates |
| `fov_um` | `(pixel_edges_x[1], pixel_edges_x[end], pixel_edges_y[1], pixel_edges_y[end])` | µm | `(xmin_um, xmax_um, ymin_um, ymax_um)`; validated `xmin < xmax`, `ymin < ymax` |

Reflection / augmentation is **generated internally** and used only for
alpha-shape polygon construction. The augmented set never appears in any
output artifact.

Coordinates: all µm. Origin and axes match `SMLMData` convention.

### 2a. Parameters & defaults

The v1 classifier is the **outer-polygon** classifier: build the polygon
from a multi-K density gate on the FOV-augmented set, then classify each
original emitter by point-in-polygon + distance to the polygon boundary.
Asymmetry-based per-emitter gates are NOT part of v1.

| Param | Default | Units | Stage | Stability |
|-------|---------|-------|-------|-----------|
| `K_LIST` | `[16, 128]` | — | multi-K density gate (augmented) | provisional |
| `RHO_K_THRESH` | 200 | µm⁻² | per-K tissue gate (intersection across K) | provisional |
| `ALPHA_NM` | 300 | nm | alpha-shape circumradius | provisional |
| `REFLECT_RADIUS_NM` | 1500 | nm | mirror band width inboard of truncated sides | provisional |
| `MEMBRANE_NM` | 100 | nm | band width for membrane class around outer polygon | provisional |
| `FOV_TRUNC_TOL_NM` | 150 | nm | truncation-detection tolerance | provisional |
| `METHOD` | `"outer_polygon"` | enum string | classifier method selector — `"outer_polygon"` (v1 default), `"grid_hybrid"` (opt-in density-grid membrane promotion), `"mask_carve"` (opt-in carve-only repair; see §4g), `"kde_valley"` (validated adaptive KDE gate; use the `kde_valley_params()` factory), or `"concave_refined"` (reserved, errors) | provisional |
| `GRID_PX_NM` | 50 | nm | grid cell size for `METHOD="grid_hybrid"`; ignored by outer-only v1 | provisional |
| `GRID_SMOOTH_NM` | 80 | nm | Gaussian smoothing σ for the density grid in `METHOD="grid_hybrid"` | provisional |
| `GRID_MASK_Q` | 0.03 | — | lower quantile floor for nonzero smoothed grid threshold in `METHOD="grid_hybrid"` | provisional |
| `GRID_MASK_PEAK_FRAC` | 0.26 | — | peak-relative smoothed grid threshold in `METHOD="grid_hybrid"` | provisional |
| `GRID_OUTER_BUFFER_NM` | 800 | nm | max distance from v1 outer polygon for interior→membrane promotion in `METHOD="grid_hybrid"` | provisional |
| `CONCAVITY_METRIC_BUFFER_NM` | 2000 | nm | buffer around the outer polygon in which `compute_concavity_metric` evaluates suspects; does not affect classification | provisional |
| `MASK_CARVE_SIGMA_UM` | 0.080 | µm | KDE Gaussian σ for the density grid used by `METHOD="mask_carve"` | provisional |
| `MASK_CARVE_K_NOISE` | 3.0 | — | multiplier on the Otsu-estimated noise floor; threshold = `K_NOISE × noise_floor` | provisional |
| `MASK_CARVE_PIXEL_UM` | 0.040 | µm | grid pixel pitch for `METHOD="mask_carve"` | provisional |
| `MASK_CARVE_MIN_COMPONENT_FRAC` | 0.05 | — | drop connected components smaller than this fraction of the largest before carving | provisional |
| `MASK_CARVE_FILL_HOLE_MAX_UM2` | 0.5 | µm² | fill internal mask holes up to this area before carving (preserves legitimate large voids) | provisional |
| `KDE_SIGMA_NM` | 150 | nm | Gaussian-KDE bandwidth σ for `METHOD="kde_valley"` (validated A431 dSTORM value) | validated |
| `KDE_RMAX_SIGMA` | 3.0 | — | KDE range-query cutoff in units of σ | validated |
| `KDE_VALLEY_NBINS` | 140 | — | log-density histogram bins for the valley threshold | validated |
| `KDE_VALLEY_FLOORFRAC` | 0.05 | — | left-base cutoff as a fraction of the cell-mode peak | validated |
| `KDE_VALLEY_SMOOTH` | 4 | bins | ±window for histogram smoothing before valley search | validated |
| `FOOTPRINT_BIN_UM` | 0.2 | µm | raster bin for the footprint fill (`METHOD="kde_valley"`) | validated |
| `FOOTPRINT_CLOSING_PX` | 3 | px | morphological closing radius to seal thin necks before hole-fill | validated |
| `ENCLOSURE_BIN_UM` | 0.2 | µm | raster bin for the 8-ray enclosure reclass | validated |
| `ENCLOSURE_MIN_HITS` | 6 | of 8 | min rays hitting cell tissue to fold a background point into `interior` | validated |

All parameter keys uppercase. Defaults will move as we tune; callers pinning
a specific parameter set should record `params_used` from `params.json`
(which records all resolved method and grid parameters as well).

---

## 3. Output layout

Nested per-cell leaf. Inside the leaf, filenames are simple and don't
repeat `<condition>` / `<cell>`.

```
<out_dir>/
  <condition>/                # e.g. RGY
    <cell>/                   # e.g. cell_01
      classified.tsv          # §4a (stable)
      polygon_loops.tsv       # §4b (stable)
      loop_diagnostics.csv    # §4c (stable cols; type provisional)
      params.json             # §4d (stable, additive only)
      manifest.json           # §4e (stable)
      classified.png          # diagnostic, written iff --renders
      loop_overlay.png        # diagnostic, written iff --renders
```

`manifest.json` is the index — consumers read it first to discover artifact
paths and schema versions, never the directory listing.

---

## 4. Output artifacts

### 4a. `classified.tsv` — STABLE

Per-emitter classification. One row per original emitter in input order.

```
# schema_version: 2
# condition: RGY
# cell: cell_01
# n_emitters: 290500
# coord_units: um
emitter_id  x_um       y_um       class      inside_outer  in_cell  dist_to_outer_um
1           4.13027    8.91204    interior   1             1        0.412
...
```

| Column | Type | Meaning |
|--------|------|---------|
| `emitter_id` | int (1-based) | row index — `emitter_id == i` joins back to `smld.emitters[i]` |
| `x_um`, `y_um` | float | echoed for joinability |
| `class` | enum string | one of `outside`, `membrane`, `interior` |
| `inside_outer` | 0/1 | **geometric** containment inside `loop_id == 1` (outer) polygon |
| `in_cell` | 0/1 | **topological** cell membership, `== (class != "outside")` (added in schema_version 2). Equals `inside_outer` for every method except `METHOD="kde_valley"`, where the enclosure stage folds enclosed background into `interior` |
| `dist_to_outer_um` | float | min perpendicular distance to outer polygon edges; `NaN` if `inside_outer == 0` |

**v1 class semantics (outer-only decision)**:

- `outside` ⇔ `inside_outer == 0`
- `membrane` ⇔ `inside_outer == 1 AND dist_to_outer_um < MEMBRANE_NM/1000`
- `interior` ⇔ `inside_outer == 1 AND dist_to_outer_um >= MEMBRANE_NM/1000`

Interior loops (`loop_id >= 2`) do **not** affect the class column in v1.
Whether they should is a v2 decision (open question in §8).

**`METHOD="kde_valley"` relaxes these biconditionals** — `class` is authoritative.
The enclosure stage reclassifies background points enclosed by the cell to
`interior` while leaving `inside_outer`/`dist_to_outer_um` strictly geometric, so
`interior ⊇ {inside_outer == 1 AND dist >= MEMBRANE}` and the enclosure-recovered
set is exactly `class == "interior" AND inside_outer == 0` (those have
`dist_to_outer_um == NaN`). `membrane` stays the band around the geometric outer
polygon only. Downstream interior filters should read `class` (or `in_cell` for
membership), never `inside_outer`.

Class invariants:

- Classes partition the input set (no nulls, no duplicates).
- Order is identical to input.

### 4b. `polygon_loops.tsv` — STABLE

All boundary loops produced by alpha-shape on the augmented set, sorted by
`abs(area)` descending. `loop_id == 1` is the outer mosaic boundary;
`loop_id >= 2` are interior / hole / reflection-space loops.

```
# schema_version: 1
# alpha_nm: 300
# reflect_radius_nm: 1500
# loop_count: 35
loop_id  vertex_id  x_um       y_um
1        1          23.70399   11.39046
...
```

Vertex coordinates **may lie outside the FOV** when they came from mirrored
emitters. Consumers wanting in-FOV-only vertices should filter against the
FOV bounds (recorded in `params.json` as `fov_um`).

### 4c. `loop_diagnostics.csv` — schema_version 2 on this branch

Per-loop summary. Includes **`loop_id == 1`** as well as all interior loops.
A `# schema_version: 2` header comment line precedes the column header.

```
# schema_version: 2
loop_id,vertex_count,area_um2,n_emitters_inside,frac_in_fov,frac_dense,median_rhoK,used_in_outer,heuristic_type
1,1060,606.49,288910,0.66,0.91,231,true,outer
2,52,4.74,314,0.50,0.44,98,false,fov_crossing
3,38,3.29,418,1.00,1.00,219,false,interior_dense
...
```

Column order is fixed: stable columns first (in the order below),
`used_in_outer` second-to-last, `heuristic_type` last.

| Column | Type | Stability | Meaning |
|--------|------|-----------|---------|
| `loop_id` | int | stable | matches `polygon_loops.tsv` |
| `vertex_count` | int | stable | number of vertices in the loop |
| `area_um2` | float | stable | `abs(polygon_area)` |
| `n_emitters_inside` | int | stable | originals strictly inside (point-in-polygon) |
| `frac_in_fov` | float in [0,1] | stable | fraction of vertices inside `fov_um` |
| `frac_dense` | float in [0,1] | stable | fraction of vertices with ρ_K(K=128) ≥ `RHO_K_THRESH`, evaluated against the **originals-only** KDTree |
| `median_rhoK` | float | stable | median ρ_K(K=128) across loop vertices (originals KDTree) |
| `used_in_outer` | bool | stable (added in schema 2) | true iff this loop participates in the `inside_outer` decision. v1 outer-polygon and `grid_hybrid`: only `loop_id == 1`. For `METHOD == "mask_carve"`, `loop_id == 1` is still marked `used_in_outer = true` because the alpha outer is the **source envelope** for the carve — but the EFFECTIVE classification boundary is the carve polygon recorded in `effective_outer.tsv`, not `loop_id == 1` of `polygon_loops.tsv`. Future methods may promote auxiliary loops. |
| `heuristic_type` | enum string | name & presence stable; values & thresholds PROVISIONAL | human-diagnostic loop label (see below) |

The `manifest.json` `artifacts.loop_diagnostics_csv.schema_version` is
bumped to `2` accordingly. Manifest top-level `schema_version` remains
`1` (additive change to one artifact).

**`heuristic_type` (column name stable, values PROVISIONAL — do not branch on it)**:

A heuristic label for human reading; thresholds may move:

- `outer` — `loop_id == 1`
- `reflection_noise` — `frac_in_fov == 0.0`
- `fov_crossing` — `0.0 < frac_in_fov < 1.0`
- `interior_dense` — `frac_in_fov == 1.0 AND frac_dense >= 0.8`
- `interior_sparse` — `frac_in_fov == 1.0 AND frac_dense < 0.8`

@genmab should consume the raw columns and threshold themselves if they
need a partition. Promotion of a loop-label ontology to stable status
would happen via a separate schema bump.

### 4d. `params.json` — STABLE schema (additive only)

Provenance contract for run comparison.

```json
{
  "schema_version": 1,
  "git_sha": "858cd25...",
  "git_status_clean": true,
  "timestamp_utc": "2026-05-04T20:43:11Z",
  "input": {
    "smld_path": "/mnt/.../Cell_01_bagol/08_bagol/smld_bagol.jld2",
    "smld_mtime_utc": "2025-06-04T...",
    "smld_size_bytes": 12345678
  },
  "n_emitters": 290500,
  "n_reflected": 33860,
  "fov_um": [0.0, 25.04, 0.0, 25.04],
  "truncated_sides": {"L": true, "R": true, "B": true, "T": true},
  "params": {
    "K_LIST": [16, 128], "RHO_K_THRESH": 200, "ALPHA_NM": 300,
    "REFLECT_RADIUS_NM": 1500, "MEMBRANE_NM": 100, "FOV_TRUNC_TOL_NM": 150,
    "METHOD": "outer_polygon",
    "GRID_PX_NM": 50, "GRID_SMOOTH_NM": 80,
    "GRID_MASK_Q": 0.03, "GRID_MASK_PEAK_FRAC": 0.26,
    "GRID_OUTER_BUFFER_NM": 800,
    "CONCAVITY_METRIC_BUFFER_NM": 2000
  },
  "runtime_s": 612.4
}
```

`input.smld_path` / `input.smld_mtime_utc` / `input.smld_size_bytes` are
omitted (or `null`) when the core API is called with raw arrays and no
SMLD path is available.

Parameter defaults are **provisional** and may move; the schema (keys,
nesting) is stable except via explicit `schema_version` bumps.

### 4e. `manifest.json` — STABLE

Single index point. Consumers read this first.

```json
{
  "schema_version": 1,
  "condition": "RGY",
  "cell": "cell_01",
  "out_dir": "/abs/path/to/out_dir",
  "leaf_dir": "/abs/path/to/out_dir/RGY/cell_01",
  "artifacts": {
    "classified_tsv":     {"path": "classified.tsv",     "schema_version": 1},
    "polygon_loops_tsv":  {"path": "polygon_loops.tsv",  "schema_version": 1},
    "loop_diagnostics_csv": {"path": "loop_diagnostics.csv", "schema_version": 2},
    "params_json":        {"path": "params.json",        "schema_version": 1},
    "classified_png":     {"path": "classified.png",     "written": false, "schema_version": null},
    "loop_overlay_png":   {"path": "loop_overlay.png",   "written": false, "schema_version": null},
    "effective_outer_tsv":         {"path": "effective_outer.tsv",         "written": false, "schema_version": null},
    "mask_carve_diagnostic_json":  {"path": "mask_carve_diagnostic.json",  "written": false, "applied": false, "schema_version": null}
  },
  "timestamp_utc": "2026-05-04T20:43:11Z"
}
```

Paths in `artifacts.*.path` are relative to `leaf_dir`. PNG entries always
appear; `written` is `false` when `--renders` was not passed. The
`effective_outer_tsv` and `mask_carve_diagnostic_json` entries always
appear (additive under manifest schema 1); `written = true` /
`schema_version = 1` only when `METHOD == "mask_carve"`. The
`applied` field on `mask_carve_diagnostic_json` mirrors the diagnostic's
own `applied` flag (`false` when carve fell back to v1).

### 4f. Renders — DIAGNOSTIC / PROVISIONAL

`classified.png` (3-class scatter) and `loop_overlay.png` (interior loops
colored by `frac_dense`) are diagnostic outputs only. Visual styling, marker
sizes, colors, axis presence are **not** part of any contract. @genmab
should consume the TSV/CSV/JSON artifacts, not the PNGs.

### 4g. `mask_carve` artifacts — PROVISIONAL / OPT-IN

These artifacts are emitted **only** when `METHOD == "mask_carve"`. The
`manifest.json` lists both entries unconditionally with `written` and
`schema_version`; for non-`mask_carve` methods `written = false` and
`schema_version = null`.

**`effective_outer.tsv`** (schema_version 1) — the polygon actually used
for `inside_outer` / `dist_to_outer` decisions. For `mask_carve` this is
the carved boundary, which differs from `polygon_loops.tsv` `loop_id == 1`
(that file always records the alpha-shape outer for provenance).

```
# schema_version: 1
# method: mask_carve
# vertex_count: <n>
vertex_id  x_um  y_um  method
1          1.234 5.678 mask_carve
…
```

**`mask_carve_diagnostic.json`** (schema_version 1) — per-call carve
diagnostic. `applied = false` indicates the carve was attempted but fell
back to the v1 alpha outer (degenerate density grid, empty intersection,
or polygonization failure); `fallback_reason` carries a human-readable
tag. Areas come from raster integration over the FOV at
`MASK_CARVE_PIXEL_UM` pitch. `carve_only_area_um2 ≈ 0` is an invariant
(carve ⊆ v1 by construction; small values are rasterization roundoff).

```json
{
  "schema_version": 1,
  "method": "mask_carve",
  "applied": true,
  "fallback_reason": "",
  "params": { "MASK_CARVE_SIGMA_UM": 0.08, "MASK_CARVE_K_NOISE": 3.0,
              "MASK_CARVE_PIXEL_UM": 0.04,
              "MASK_CARVE_MIN_COMPONENT_FRAC": 0.05,
              "MASK_CARVE_FILL_HOLE_MAX_UM2": 0.5 },
  "v1_polygon_area_um2": 250.5,
  "carve_polygon_area_um2": 244.8,
  "area_delta_um2": -5.7,
  "v1_only_area_um2": 5.7,
  "carve_only_area_um2": 0.0,
  "med_v1_carve_distance_um": 1.21,
  "p95_v1_carve_distance_um": 1.71,
  "n_holes_filled": 12,
  "n_holes_preserved": 0,
  "n_carve_polygon_pts": 2345
}
```

**Carve-only limitation.** `mask_carve` only carves the v1 outer polygon
inward; it cannot recover membrane that lies outside v1. Regions where
v1 already underestimates the cell are unchanged. Provisional / opt-in
on this branch — defaults are taken from synthetic stop-condition gates
in `dev/scripts/mask_contour_v3.jl` and `mask_contour_v31.jl` (no
real-cell tuning); behavior on real cells is documented in the v3.1 dev
report.

---

## 5. Class invariants (summary)

- `class` in `classified.tsv` is the canonical per-emitter answer.
- `outside ∪ membrane ∪ interior == all input emitters`, intersections empty.
- `METHOD="outer_polygon"` decision: outer loop only. Interior loops are diagnostic.
- `METHOD="grid_hybrid"` preserves the outer-loop `outside`/`interior`
  topology, then promotes only v1 `interior` emitters to `membrane` when
  they lie both on the local density-grid boundary and within
  `GRID_OUTER_BUFFER_NM` of the v1 outer polygon. It never demotes
  `membrane`, never changes `outside`, and does not use diagnostic
  interior loops for class decisions.
- `METHOD="mask_carve"` replaces the **effective** outer polygon used for
  `inside_outer` / `dist_to_outer` with a carved subset of v1 (alpha
  outer remains in `loops[1]` and `polygon_loops.tsv` for provenance).
  Carve ⊆ v1 by construction. May reclassify v1 `interior` / `membrane`
  emitters in carved-away regions to `outside`, and may demote
  `membrane` to `interior` near a carve boundary; never adds emitters
  outside v1 to the cell.
- `METHOD="kde_valley"` (validated genmab dSTORM gate) replaces the **density
  gate** upstream of the polygon: a continuous Gaussian-KDE density is thresholded
  at the background/cell valley (per-FOV adaptive — handles the ~6× MAP-N density
  spread with no per-cell tuning), footprint-filled, then the v1 outer-polygon
  geometry runs on the footprint subset, and an enclosure pass folds background
  points enclosed by the cell into `interior`. `class` is authoritative and
  includes the enclosure-recovered interiors; `inside_outer`/`dist_to_outer_um`
  stay strictly geometric (the enclosure-recovered set is `class == "interior" AND
  inside_outer == 0`); `in_cell == (class != "outside")` carries topological
  membership. Use the `kde_valley_params()` factory (validated defaults σ=150 nm,
  α=600 nm, reflect=1500 nm, membrane=100 nm) — the struct default `ALPHA_NM` is
  300, so a raw constructor would under-alpha. dSTORM (A431/HeLa) path only.

---

## 6. Provenance / run comparison

`params.json` is the single source of truth for "how was this run
configured?". @genmab should diff this file across runs.

`manifest.json` is the single source of truth for "what artifacts exist
and where?". Don't list directories.

---

## 7. Stable vs provisional vs internal

| Stable (won't break w/o schema bump) | Provisional (may move) | Internal (do not depend on) |
|--------------------------------------|------------------------|-----------------------------|
| Filenames in §3 leaf (`classified.tsv`, `polygon_loops.tsv`, `loop_diagnostics.csv`, `params.json`, `manifest.json`) | Default parameter values | Helper scripts under `dev/scripts/` (`polygon_reflected.jl`, `tissue_feature_histograms.jl`, `polygon_reflected_diagnostics.jl`, etc.) |
| `classified.tsv` columns & class labels | `loop_diagnostics.type` heuristic | Render visual styling |
| `polygon_loops.tsv` schema | Inclusion of optional diagnostic columns beyond the core | Augmented-set / reflection internals |
| `loop_diagnostics.csv` core columns + `heuristic_type` column name & presence | Render PNG existence / contents; `heuristic_type` values & thresholds | Module names below `EdgeClassify` public surface |
| `params.json` keys (additive only) | v1 class decision rule (outer-only) — may extend in v2 |  |
| `manifest.json` schema (additive only) |  |  |
| Class invariants in §5 |  |  |

Schema bumps to stable artifacts: bump `schema_version` and document in
this file; never break old consumers without explicit notice.

---

## 8. Resolved decisions (signed off by @codex-cluster)

1. **v2 interior-loop class contribution: deferred**. v1 `class` is
   outer-only. Any future internal void/membrane semantics ship as a
   schema bump or new columns/artifacts after biological review — not as
   a silent change.
2. **Heuristic loop label: keep, renamed `heuristic_type`**, placed last
   in `loop_diagnostics.csv`. Column name and presence are stable; values
   and thresholds are provisional. Consumers branch on raw metrics only.
3. **`fov_um` order: pinned `(xmin_um, xmax_um, ymin_um, ymax_um)`**.
   Used everywhere (API docs, `params.json`, validation). Implementation
   validates `xmin < xmax` and `ymin < ymax`.
4. **`params.toml` schema**: same keys as `params.json.params`, uppercase
   to match existing scripts. Unknown keys → error. Missing keys → defaults.
   `params.json` records the **fully resolved** parameter set used.

---

## 9. Implementation target

After this doc is locked, implement `SMLMClustering.EdgeClassify` with:

- coordinate/FOV core function
- SMLD adapter
- CLI script (`dev/scripts/edge_classify.jl`)
- artifacts + manifest writers per §3 / §4
- tests / smoke checks for schema and class invariants
- renders optional behind `--renders`

---

## Corrections from v0.1 (acknowledged)

- **Artifact paths**: previously documented as `output/...` relative to
  whatever cwd the script ran in. v1 makes `--out` authoritative; absolute
  paths only.
- **Population count**: 20 interior_dense + 9 fov_crossing + 5 reflection_noise
  = 34 interior loops (not 21 + 9 + 5). Run 1 had 20 interior_dense; the
  +14 in Run 2 = 9 fov_crossing + 5 reflection_noise, all attributable to
  the reflection step.
