# Local-contrast density feature for the cluster_statistics interface.
#
# Per-emitter feature equal to the kNN log-density at a fine scale minus the
# median kNN log-density over a larger neighborhood. This captures *local*
# elevation in density relative to nearby baseline, which is robust to global
# density gradients across a cell that confound absolute-density thresholds.
#
# Mechanism: the absolute-density classifier (Otsu / GMM on log ρ_k) calls a
# point "structure" whenever its local density exceeds a single global cutoff.
# When baseline density varies across the field of view (cell-edge thinning,
# illumination falloff, biological gradient) the cutoff is wrong on at least
# one side. Subtracting a coarse local baseline cancels the gradient: the
# feature only fires when the point is denser than its surroundings, not just
# denser than the cell-wide average.

import Statistics
import NearestNeighbors
import Base.Threads: @threads, threadid, nthreads, maxthreadid

"""
    LocalContrastFeature(; density_k=200, background_k=2000,
                          use_3d=false, per_dataset=false)

Per-emitter local-density-contrast feature.

# Algorithm
1. Build a `KDTree` over the (per-dataset or pooled) emitter coordinates.
2. For each emitter `i`, compute the kNN log-density
   `f_i = log(density_k / (π · r_k²))` where `r_k` is the distance to the
   `density_k`-th nearest neighbor (excluding self). The kNN ball area
   normalization gives a density estimate in coordinate-units⁻²; the log keeps
   the downstream median/threshold comparisons in a comparable scale across a
   cell.
3. For each emitter `i`, compute the median of `f_j` over its `background_k`
   nearest neighbors (excluding self) — the *local baseline* at a coarser
   spatial scale. Median is used (not mean) because it tolerates a small
   number of locally-elevated neighbors without dragging the baseline up.
4. Contrast `c_i = f_i − median(f_j over j in k_bg-NN of i)`. Positive contrast
   means point `i` is denser than its local surroundings.

# Fields
- `density_k::Int = 200`: fine-scale neighborhood for the per-point log-density.
  Sets the spatial scale of the *signal*; for SMLM cells, ~80–200 nm typically.
- `background_k::Int = 2000`: coarse-scale neighborhood for the local baseline
  (must be > `density_k`). Sets the spatial scale of the *baseline*; should be
  large enough to span the structures you want to detect (typically 1–2 μm).
- `use_3d::Bool = false`: if `true`, build the KDTree over `(x, y, z)`.
- `per_dataset::Bool = false`: if `true`, compute per dataset independently
  (each dataset gets its own KDTree). Per-emitter outputs stitch back into
  original emitter order in either case.

# Returned info
- `statistic = median(non-NaN contrast)` — single scalar summary. NaN when no
  group has enough points to compute the feature.
- `statistic_name = :median_local_contrast`
- `algorithm = :local_contrast`
- `extras[:contrast_per_emitter]` — `Vector{Float64}` of length `n_locs_in`,
  in original emitter order. Units: nat (natural-log scale). The caller
  thresholds this directly.
- `extras[:log_density_per_emitter]` — `Vector{Float64}` of length `n_locs_in`,
  the fine kNN log-density `f_i`. Useful for absolute-density gates that
  complement the local-contrast gate (e.g. require `f_i > q35` floor).

# Edge case handling
- Groups with `n ≤ density_k`: those emitters receive `NaN` for both contrast
  and log-density. Avoids a degenerate kNN query against a too-small set.
- `background_k ≥ n` in a group: clamped to `n − 1` for that group; the
  feature becomes "log-density minus group median," which is still
  well-defined.
- `background_k ≤ density_k`: raises `ArgumentError` at config use; the
  baseline must be coarser than the signal.
- Coincident coordinates (`r_k = 0`): the affected emitter receives `NaN`
  log-density and `NaN` contrast.

# Composition
This feature is the missing primitive for hysteresis seed-and-grow on point
clouds with non-stationary baseline density:
```julia
(_, info) = cluster_statistics(smld, LocalContrastFeature())
contrast = info.extras[:contrast_per_emitter]
fine = info.extras[:log_density_per_emitter]
fine_floor = quantile(filter(isfinite, fine), 0.35)
seed = isfinite.(contrast) .& (contrast .> 0.25) .& (fine .> fine_floor)
support = isfinite.(contrast) .& (contrast .> -0.05) .& (fine .> fine_floor)
(smld_out, _) = cluster(smld, PointHysteresisConfig(graph_k=12, min_points=150);
                        seed=seed, support=support)
```

See also: [`AbstractStatisticsConfig`](@ref), [`ClusterStatisticsInfo`](@ref).
"""
Base.@kwdef struct LocalContrastFeature <: AbstractStatisticsConfig
    density_k::Int = 200
    background_k::Int = 2000
    use_3d::Bool = false
    per_dataset::Bool = false
end

function _local_contrast_group!(contrast::Vector{Float64},
                                log_density::Vector{Float64},
                                emitters,
                                idxs::AbstractVector{Int},
                                density_k::Int,
                                background_k::Int,
                                use_3d::Bool)
    n = length(idxs)
    if n <= density_k
        # Cannot compute the feature — leave NaN entries for these emitters.
        return
    end
    sub = view(emitters, idxs)
    X = _coords_matrix(sub, use_3d)
    tree = NearestNeighbors.KDTree(X)
    fine_local = Vector{Float64}(undef, n)

    @threads for j in 1:n
        _, dists = NearestNeighbors.knn(tree, view(X, :, j), density_k + 1, true)
        rk = dists[end]
        fine_local[j] = rk > 0 ? log(density_k / (π * rk^2)) : NaN
    end

    k_bg = min(background_k, n - 1)
    # Index buffers by `threadid()`, which on Julia 1.10+ can exceed
    # `nthreads()` because of interactive / default thread pools — use
    # `maxthreadid()` to size the buffer table safely.
    nbuf = max(1, maxthreadid())
    bufs = [Vector{Float64}(undef, k_bg) for _ in 1:nbuf]

    @threads for j in 1:n
        nbrs, _ = NearestNeighbors.knn(tree, view(X, :, j), k_bg + 1, true)
        buf = bufs[threadid()]
        m = 0
        @inbounds for nb in nbrs
            nb == j && continue
            v = fine_local[nb]
            isfinite(v) || continue
            m += 1
            buf[m] = v
        end
        global_i = idxs[j]
        log_density[global_i] = fine_local[j]
        if m == 0 || !isfinite(fine_local[j])
            contrast[global_i] = NaN
        else
            contrast[global_i] = fine_local[j] - Statistics.median!(@view buf[1:m])
        end
    end
    return
end

function cluster_statistics(smld::SMLMData.BasicSMLD, cfg::LocalContrastFeature)
    t0 = time_ns()
    n_in = length(smld.emitters)

    cfg.background_k > cfg.density_k || throw(ArgumentError(
        "LocalContrastFeature requires background_k > density_k " *
        "(got density_k=$(cfg.density_k), background_k=$(cfg.background_k)). " *
        "The local baseline must be coarser than the signal."))
    cfg.density_k >= 1 || throw(ArgumentError(
        "LocalContrastFeature density_k must be ≥ 1 (got $(cfg.density_k))."))

    contrast = fill(NaN, n_in)
    log_density = fill(NaN, n_in)

    groups = _group_by_dataset(smld, cfg.per_dataset)
    for idxs in groups
        _local_contrast_group!(contrast, log_density, smld.emitters, idxs,
                               cfg.density_k, cfg.background_k, cfg.use_3d)
    end

    valid = filter(isfinite, contrast)
    median_contrast = isempty(valid) ? NaN : Statistics.median(valid)

    extras = Dict{Symbol,Any}(
        :contrast_per_emitter => contrast,
        :log_density_per_emitter => log_density,
    )

    info = ClusterStatisticsInfo(
        n_in,
        median_contrast,
        :median_local_contrast,
        :local_contrast,
        (time_ns() - t0) / 1e9,
        extras,
    )
    return smld, info
end
