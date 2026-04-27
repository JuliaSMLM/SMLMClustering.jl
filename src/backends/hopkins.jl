# Hopkins-statistic backend for the cluster_statistics interface.
#
# Computes the Hopkins statistic H, a measure of clustering tendency on an
# unlabeled point set. Definition:
#
#   For a sample of m reference points uniformly drawn from the data's
#   bounding box, let u_i be each reference point's NN distance to the data.
#   For a sample of m real data points (drawn without replacement from the
#   data set), let w_i be each sampled point's NN distance to the OTHER
#   real data points (excluding itself).
#
#   H = sum(u_i^d) / (sum(u_i^d) + sum(w_i^d))    (d = 2 or 3)
#
#   - H ≈ 0.5 ⇒ the data is statistically indistinguishable from uniform
#   - H → 1.0 ⇒ strong clustering tendency
#   - H → 0.0 ⇒ regular / lattice-like
#
# Per-dataset computation reports the mean H across datasets in `statistic`
# and the full per-dataset vector in `extras[:hopkins_per_dataset]`.

"""
    HopkinsConfig(; n_samples=20, random_repeats=1, seed=nothing,
                    use_3d=false, per_dataset=true)

Configuration for the Hopkins-statistic spatial-randomness test.

# Algorithm
For each repeat:
1. Draw `n_samples` reference points uniformly from the per-group bounding
   box; let `u_i` = NN distance from each reference point to the data.
2. Draw `n_samples` real points without replacement from the data; let
   `w_i` = NN distance from each sampled point to the OTHER real points
   (excluding itself).
3. With `d` = 2 (or 3 if `use_3d=true`), `H = Σuᵢ^d / (Σuᵢ^d + Σwᵢ^d)`.

Repeats are averaged.

# Fields
- `n_samples::Int = 20`: number of reference / sampled points per repeat.
  Must satisfy `n_samples ≤ n_points` per group.
- `random_repeats::Int = 1`: number of independent repeats to average. Higher
  values reduce variance at linear cost.
- `seed::Union{Int,Nothing} = nothing`: when set, seeds an internal `Xoshiro`
  for reproducibility. When `nothing`, uses the global RNG.
- `use_3d::Bool = false`: include the z-coordinate (and use `d=3` in the formula).
- `per_dataset::Bool = true`: when `true`, compute Hopkins per dataset; the
  reported `statistic` is the mean across datasets and the full per-dataset
  vector is placed in `extras[:hopkins_per_dataset]`. When `false`, all
  emitters are pooled and a single H is returned.

# Interpretation
- `H ≈ 0.5`: data is consistent with uniform spatial randomness (Poisson)
- `H → 1.0`: strong clustering tendency
- `H → 0.0`: anti-clustering / regular spacing

# Example
```julia
cfg = HopkinsConfig(n_samples = 50, random_repeats = 5, seed = 1)
(_, info) = cluster_statistics(smld, cfg)
println("Hopkins H = \$(round(info.statistic, digits=3))")
```

See also: [`AbstractStatisticsConfig`](@ref), [`ClusterStatisticsInfo`](@ref),
[`cluster_statistics`](@ref).
"""
Base.@kwdef struct HopkinsConfig <: AbstractStatisticsConfig
    n_samples::Int = 20
    random_repeats::Int = 1
    seed::Union{Int,Nothing} = nothing
    use_3d::Bool = false
    per_dataset::Bool = true
end

# Compute Hopkins H on a single d×n coordinate matrix using a supplied RNG.
# Returns NaN when `n < 2` or `n_samples > n` (caller decides how to aggregate).
function _hopkins_one_group(X::Matrix{Float64}, n_samples::Int, repeats::Int, rng)
    d, n = size(X)
    (n >= 2 && n_samples <= n) || return NaN

    # Bounding box of the group.
    lo = Vector{Float64}(undef, d)
    hi = Vector{Float64}(undef, d)
    @inbounds for k in 1:d
        col = @view X[k, :]
        lo[k] = minimum(col)
        hi[k] = maximum(col)
    end
    extent = hi .- lo
    # Degenerate bbox (all points coincident in some axis): zero-volume box,
    # uniform sampling is degenerate. Return NaN — callers can interpret.
    any(==(0.0), extent) && return NaN

    tree = KDTree(X)
    Hs = Vector{Float64}(undef, repeats)

    for r in 1:repeats
        u_sum = 0.0
        w_sum = 0.0

        # 1. Reference points: uniform in bbox; NN to data.
        for _ in 1:n_samples
            ref = Vector{Float64}(undef, d)
            @inbounds for k in 1:d
                ref[k] = lo[k] + extent[k] * rand(rng)
            end
            _, dists = NearestNeighbors.knn(tree, ref, 1, true)
            u_sum += dists[1]^d
        end

        # 2. Sampled real points: drawn without replacement; NN to OTHER real points.
        sample_idx = randperm(rng, n)[1:n_samples]
        for j in sample_idx
            pt = @view X[:, j]
            # k=2 because the nearest neighbor of a data point in its own tree
            # is itself (distance 0); we want the second-nearest.
            _, dists = NearestNeighbors.knn(tree, Vector(pt), 2, true)
            w_sum += dists[2]^d
        end

        denom = u_sum + w_sum
        Hs[r] = denom == 0.0 ? NaN : u_sum / denom
    end

    return sum(Hs) / repeats
end

function cluster_statistics(smld::SMLMData.BasicSMLD, cfg::HopkinsConfig)
    t0 = time_ns()
    n_in = length(smld.emitters)
    cfg.n_samples >= 1 ||
        throw(ArgumentError("HopkinsConfig.n_samples must be ≥ 1 (got $(cfg.n_samples))"))
    cfg.random_repeats >= 1 ||
        throw(ArgumentError("HopkinsConfig.random_repeats must be ≥ 1 (got $(cfg.random_repeats))"))

    rng = cfg.seed === nothing ? Random.default_rng() : Xoshiro(cfg.seed)

    groups = _group_by_dataset(smld, cfg.per_dataset)

    extras = Dict{Symbol,Any}()
    if cfg.per_dataset
        per_ds = Vector{Float64}(undef, length(groups))
        for (gi, idxs) in pairs(groups)
            if isempty(idxs)
                per_ds[gi] = NaN
                continue
            end
            sub = view(smld.emitters, idxs)
            X = _coords_matrix(sub, cfg.use_3d)
            per_ds[gi] = _hopkins_one_group(X, cfg.n_samples, cfg.random_repeats, rng)
        end
        extras[:hopkins_per_dataset] = per_ds
        # Mean across non-NaN datasets; if all NaN, statistic is NaN.
        valid = filter(!isnan, per_ds)
        H = isempty(valid) ? NaN : sum(valid) / length(valid)
    else
        # Single pooled computation across all emitters.
        if isempty(smld.emitters)
            H = NaN
        else
            X = _coords_matrix(smld.emitters, cfg.use_3d)
            H = _hopkins_one_group(X, cfg.n_samples, cfg.random_repeats, rng)
        end
    end

    info = ClusterStatisticsInfo(
        n_in,
        H,
        :hopkins,
        :hopkins,
        (time_ns() - t0) / 1e9,
        extras,
    )
    # Pass-through SMLD reference per V10 (no copy, no mutation).
    return smld, info
end
