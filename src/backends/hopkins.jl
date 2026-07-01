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
- `region = nothing`: observation window for the uniform **reference** points
  (2D only). Hopkins is window-sensitive — sampling references over the data
  bounding box makes data that is uniform inside a non-convex boundary read as
  *falsely* clustered. Options: `nothing` → the data bounding box (default); a
  polygon `Vector{NTuple{2,Float64}}` → references are rejection-sampled inside it;
  `:metadata` → use the polygon at `smld.metadata["edge_outer_polygon"]` (written
  by `classify_emitters` — the pipeline channel); `Dict(dataset_id => polygon)` →
  one polygon per dataset (`per_dataset = true`). Incompatible with `use_3d = true`.

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
    region::Union{Nothing,Symbol,Vector{NTuple{2,Float64}},Vector{CellPolygon},Dict{Int,Vector{NTuple{2,Float64}}}} = nothing
end

# Reference-region membership + bounding box, for either a single polygon or a
# multi-cell mask used as a Hopkins observation window.
_region_contains(region::Vector{NTuple{2,Float64}}, x, y) = _point_in_polygon(x, y, region)
_region_contains(region::AbstractVector{CellPolygon}, x, y) = in_region(x, y, region)

function _region_bbox(region::Vector{NTuple{2,Float64}})
    lox = loy = Inf; hix = hiy = -Inf
    @inbounds for (vx, vy) in region
        lox = min(lox, vx); hix = max(hix, vx); loy = min(loy, vy); hiy = max(hiy, vy)
    end
    return lox, hix, loy, hiy
end
function _region_bbox(region::AbstractVector{CellPolygon})
    lox = loy = Inf; hix = hiy = -Inf
    @inbounds for cell in region
        for (vx, vy) in cell.outer
            lox = min(lox, vx); hix = max(hix, vx); loy = min(loy, vy); hiy = max(hiy, vy)
        end
    end
    return lox, hix, loy, hiy
end

# Compute Hopkins H on a single d×n coordinate matrix using a supplied RNG.
# `region`, when given, is a 2D polygon: reference points are rejection-sampled
# inside it (the correct observation window) instead of the data bounding box, so
# data that is uniform but confined to a non-convex domain no longer reads as
# falsely clustered. Returns NaN when `n < 2`, `n_samples > n`, the sampling
# envelope is degenerate, or (region mode) a reference point cannot be placed
# inside the polygon within the attempt cap.
function _hopkins_one_group(X::Matrix{Float64}, n_samples::Int, repeats::Int, rng;
                            region::Union{Nothing,Vector{NTuple{2,Float64}},
                                          Vector{CellPolygon}} = nothing)
    d, n = size(X)
    (n >= 2 && n_samples <= n) || return NaN

    # Reference-sampling envelope: the data bounding box by default, or the region's
    # bounding box (with rejection back into the region) in region mode. The region
    # is a single polygon or a multi-cell mask; 2D — use_3d is guarded upstream, so
    # d == 2 whenever region !== nothing.
    use_region = region !== nothing
    lo = Vector{Float64}(undef, d)
    hi = Vector{Float64}(undef, d)
    if use_region
        lox, hix, loy, hiy = _region_bbox(region)
        lo[1] = lox; lo[2] = loy; hi[1] = hix; hi[2] = hiy
    else
        @inbounds for k in 1:d
            col = @view X[k, :]
            lo[k] = minimum(col)
            hi[k] = maximum(col)
        end
    end
    extent = hi .- lo
    # Degenerate envelope (zero-extent bbox, or zero-area polygon bbox): uniform
    # sampling is undefined. Return NaN — callers can interpret.
    any(==(0.0), extent) && return NaN

    # Cap rejection attempts per reference point so a thin/near-empty polygon can't
    # hang; on exhaustion the group's H is NaN.
    max_attempts = use_region ? 1000 : 1

    tree = KDTree(X)
    Hs = Vector{Float64}(undef, repeats)

    # Workspace for sampling n_samples indices without replacement via partial
    # Fisher-Yates. Allocated once per group.
    workspace = collect(1:n)

    for r in 1:repeats
        u_sum = 0.0
        w_sum = 0.0

        # 1. Reference points: uniform in the envelope; NN to data. In region mode,
        # reject samples that fall outside the polygon.
        for _ in 1:n_samples
            ref = Vector{Float64}(undef, d)
            if use_region
                placed = false
                for _ in 1:max_attempts
                    @inbounds ref[1] = lo[1] + extent[1] * rand(rng)
                    @inbounds ref[2] = lo[2] + extent[2] * rand(rng)
                    if _region_contains(region, ref[1], ref[2])
                        placed = true
                        break
                    end
                end
                placed || return NaN
            else
                @inbounds for k in 1:d
                    ref[k] = lo[k] + extent[k] * rand(rng)
                end
            end
            _, dists = NearestNeighbors.knn(tree, ref, 1, true)
            u_sum += dists[1]^d
        end

        # 2. Sampled real points: drawn without replacement via partial
        # Fisher-Yates on `workspace`. After the inner loop, workspace[1:n_samples]
        # holds n_samples distinct indices from 1..n.
        @inbounds for k in 1:n_samples
            j = rand(rng, k:n)
            workspace[k], workspace[j] = workspace[j], workspace[k]
        end
        @inbounds for k in 1:n_samples
            j = workspace[k]
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

# Validate a Hopkins region polygon.
function _check_hopkins_polygon(p::Vector{NTuple{2,Float64}})
    length(p) >= 3 ||
        throw(ArgumentError("HopkinsConfig.region polygon needs ≥ 3 vertices (got $(length(p)))"))
    all(t -> isfinite(t[1]) && isfinite(t[2]), p) ||
        throw(ArgumentError("HopkinsConfig.region polygon has non-finite vertices"))
    return nothing
end

# Resolve cfg.region into one polygon (or nothing) per group, aligned with
# `groups` (which _group_by_dataset returns in ascending dataset-id order).
function _resolve_hopkins_regions(smld::SMLMData.BasicSMLD, cfg::HopkinsConfig, groups)
    ng = length(groups)
    T = Union{Nothing,Vector{NTuple{2,Float64}},Vector{CellPolygon}}
    r = cfg.region
    r === nothing && return T[nothing for _ in 1:ng]
    if r === :metadata
        # Prefer the published multi-cell mask; fall back to the single dominant-cell polygon.
        cells = get(smld.metadata, "edge_cells", nothing)
        if cells isa Vector{CellPolygon} && !isempty(cells)
            return T[cells for _ in 1:ng]
        end
        poly = get(smld.metadata, "edge_outer_polygon", nothing)
        poly isa Vector{NTuple{2,Float64}} ||
            throw(ArgumentError("HopkinsConfig.region=:metadata requires " *
                "smld.metadata[\"edge_cells\"]::Vector{CellPolygon} or " *
                "smld.metadata[\"edge_outer_polygon\"]::Vector{NTuple{2,Float64}} " *
                "(run classify_emitters upstream); got " *
                (poly === nothing ? "nothing" : string(typeof(poly)))))
        _check_hopkins_polygon(poly)
        return T[poly for _ in 1:ng]
    elseif r isa Vector{CellPolygon}
        isempty(r) && throw(ArgumentError("HopkinsConfig.region MultiCellMask is empty"))
        return T[r for _ in 1:ng]
    elseif r isa Vector{NTuple{2,Float64}}
        _check_hopkins_polygon(r)
        return T[r for _ in 1:ng]
    elseif r isa Dict
        cfg.per_dataset ||
            throw(ArgumentError("HopkinsConfig.region as a Dict requires per_dataset=true"))
        ids = sort!(unique(e.dataset for e in smld.emitters))
        length(ids) == ng ||
            error("HopkinsConfig.region Dict: dataset-id count ($(length(ids))) ≠ group count ($ng)")
        out = T[nothing for _ in 1:ng]
        @inbounds for (gi, did) in enumerate(ids)
            haskey(r, did) ||
                throw(ArgumentError("HopkinsConfig.region Dict has no polygon for dataset $did"))
            _check_hopkins_polygon(r[did])
            out[gi] = r[did]
        end
        return out
    else
        throw(ArgumentError("HopkinsConfig.region: unsupported value of type $(typeof(r)); use " *
            "nothing, :metadata, a Vector{NTuple{2,Float64}}, or a Dict{Int,Vector{NTuple{2,Float64}}}"))
    end
end

function cluster_statistics(smld::SMLMData.BasicSMLD, cfg::HopkinsConfig)
    t0 = time_ns()
    n_in = length(smld.emitters)
    cfg.n_samples >= 1 ||
        throw(ArgumentError("HopkinsConfig.n_samples must be ≥ 1 (got $(cfg.n_samples))"))
    cfg.random_repeats >= 1 ||
        throw(ArgumentError("HopkinsConfig.random_repeats must be ≥ 1 (got $(cfg.random_repeats))"))
    (cfg.use_3d && cfg.region !== nothing) &&
        throw(ArgumentError("HopkinsConfig: `region` (a 2D polygon) is not supported with use_3d=true"))

    rng = cfg.seed === nothing ? Random.default_rng() : Xoshiro(cfg.seed)

    groups = _group_by_dataset(smld, cfg.per_dataset)
    regions = _resolve_hopkins_regions(smld, cfg, groups)

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
            per_ds[gi] = _hopkins_one_group(X, cfg.n_samples, cfg.random_repeats, rng;
                                            region = regions[gi])
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
            H = _hopkins_one_group(X, cfg.n_samples, cfg.random_repeats, rng;
                                   region = regions[1])
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
