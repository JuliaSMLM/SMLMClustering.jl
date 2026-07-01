# Voronoi-density backend for the cluster_statistics interface.
#
# Computes the per-emitter Voronoi cell area and corresponding local density
# ρᵢ = 1/Aᵢ on a `BasicSMLD`. Returns the per-emitter density and area
# vectors in `extras`, and the median density (across all emitters that
# received a valid Voronoi cell) as the summary `statistic`.
#
# This is the read-only sibling to the `VoronoiConfig` clustering backend:
# `VoronoiConfig` uses the same per-emitter areas to threshold "dense"
# emitters and form clusters; `VoronoiDensityConfig` exposes the underlying
# density measure directly so downstream callers (e.g. cell-structure
# masking via Otsu / GMM on log ρ) can run their own thresholding.
#
# 2D only — DelaunayTriangulation.jl does not provide 3D Voronoi (V7).

import Statistics

"""
    VoronoiDensityConfig(; use_3d=false, per_dataset=true)

Configuration for the per-emitter Voronoi-density spatial-statistic backend.

# Algorithm
1. For each group (dataset when `per_dataset=true`, all emitters otherwise),
   build the Voronoi tessellation of the emitter coordinates clipped to the
   convex hull (so every generator has a finite cell area).
2. For each emitter `i`, compute `Aᵢ` = its Voronoi cell area (μm²) and
   `ρᵢ = 1/Aᵢ` (μm⁻²).
3. Stitch the per-group vectors back into original emitter order so
   `extras[:density_per_emitter][i]` corresponds to `smld.emitters[i]`.

# Fields
- `use_3d::Bool = false`: must be `false`. 3D Voronoi tessellation is not
  supported (DelaunayTriangulation.jl is 2D only); passing `true` raises
  `ArgumentError`.
- `per_dataset::Bool = true`: when `true`, tessellate each dataset
  independently. When `false`, all emitters are pooled into one
  tessellation. Per-emitter outputs remain flat in original emitter order
  in either case.

# Returned info
- `statistic = median(non-NaN densities)` — single number summarizing
  overall density. NaN if no group produced any valid density.
- `statistic_name = :median_density`
- `algorithm = :voronoi_density`
- `extras[:density_per_emitter]` — `Vector{Float64}` of length `n_locs_in`,
  in original emitter order. Emitters in groups smaller than 3 (untessellatable)
  receive `NaN`. Units: μm⁻².
- `extras[:area_per_emitter]` — `Vector{Float64}` of length `n_locs_in`,
  in original emitter order. Same `NaN` semantics. Units: μm².

# Edge case handling (per V10 NaN-vs-throw rule)
- Groups with fewer than 3 points: those emitters receive `NaN` density and
  `NaN` area. Other groups proceed normally.
- Empty SMLD: empty per-emitter vectors, `statistic = NaN`.
- Group with exact-duplicate `(x, y)` coordinates: `ArgumentError` raised
  before triangulation (mirrors `VoronoiConfig`'s guard — duplicate
  coordinates are a boundary-input issue, not a data-shape edge case).

# Example
```julia
cfg = VoronoiDensityConfig()
(_, info) = cluster_statistics(smld, cfg)
ρ = info.extras[:density_per_emitter]   # Vector{Float64}, length == n_locs_in
A = info.extras[:area_per_emitter]
println("median density = \$(round(info.statistic, digits=2)) μm⁻²")
```

See also: [`AbstractStatisticsConfig`](@ref), [`ClusterStatisticsInfo`](@ref),
[`cluster_statistics`](@ref), [`VoronoiConfig`](@ref) (the clustering sibling).
"""
Base.@kwdef struct VoronoiDensityConfig <: AbstractStatisticsConfig
    use_3d::Bool = false
    per_dataset::Bool = true
end

function cluster_statistics(smld::SMLMData.BasicSMLD, cfg::VoronoiDensityConfig)
    t0 = time_ns()
    n_in = length(smld.emitters)
    cfg.use_3d &&
        throw(ArgumentError(
            "VoronoiDensityConfig does not support use_3d=true. " *
            "Voronoi tessellation is 2D only — see VoronoiConfig docstring " *
            "(KB V7) for the rationale."))

    # Per-emitter outputs are flat in original emitter order. Pre-fill with
    # NaN so untessellatable groups (< 3 points) end up correctly NaN'd
    # without any extra bookkeeping.
    density_per_emitter = fill(NaN, n_in)
    area_per_emitter = fill(NaN, n_in)

    groups = _group_by_dataset(smld, cfg.per_dataset)

    for idxs in groups
        n = length(idxs)
        n < 3 && continue  # those emitters keep NaN

        sub = view(smld.emitters, idxs)
        # `_voronoi_areas` raises ArgumentError on exact-duplicate (x,y) pairs
        # (mirrors voronoi.jl's guard); returns a Vector{Float64} of length n.
        areas, _ = _voronoi_areas(sub)

        @inbounds for j in 1:n
            a = areas[j]
            i = idxs[j]
            area_per_emitter[i] = a
            density_per_emitter[i] = (isfinite(a) && a > 0) ? 1.0 / a : NaN
        end
    end

    # Summary scalar: median over all emitters that received a valid Voronoi
    # cell (i.e. across non-NaN entries). If nothing is valid, NaN.
    valid_density = filter(!isnan, density_per_emitter)
    median_density = isempty(valid_density) ? NaN : Statistics.median(valid_density)

    extras = Dict{Symbol,Any}(
        :density_per_emitter => density_per_emitter,
        :area_per_emitter => area_per_emitter,
    )

    info = ClusterStatisticsInfo(
        n_in,
        median_density,
        :median_density,
        :voronoi_density,
        (time_ns() - t0) / 1e9,
        extras,
    )
    # Pass-through SMLD reference per V10 (no copy, no mutation).
    return smld, info
end
