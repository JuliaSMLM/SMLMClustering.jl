# MRF density-regime clustering backend.
#
# Adaptive-density clustering pipeline that auto-handles datasets with
# multiple density regimes (e.g. dSTORM data containing tight ~25 nm
# aggregates alongside μm-scale extended structure) without per-dataset
# parameter tuning.
#
# Pipeline (per group when `per_dataset=true`):
#   1. Per-emitter Voronoi density ρᵢ = 1/Aᵢ → log ρᵢ.
#   2. Regime assignment via:
#       (a) explicit `regime_thresholds` (binning), or
#       (b) `n_regimes`-component 1D Gaussian mixture EM on log ρ
#           → sorted ascending by mean (regime 1 = lowest density).
#      Per-emitter unary cost matrix `U[i, k] = -log(w_k * N(log_rho[i] | μ_k, σ_k²))`.
#   3. Multi-class Potts MRF refinement via Iterated Conditional Modes (ICM)
#      over the Delaunay neighbor graph (default) or a k-NN graph.
#      Auto-tuned smoothness λ when not user-overridden.
#   4. Connected components on the foreground (regime ≥ 2) over the same
#      neighbor graph, with `min_points` size filter.
#
# Convention: regime 1 (lowest density) is treated as noise/background;
# regimes 2..n_regimes form the foreground from which clusters are
# extracted via connected components on the neighbor graph.
#
# Output:
# - `emitter.id`: 0 (noise) or cluster id 1..K (per group, V3 namespacing).
# - `metadata["mrf_regime_per_emitter"]::Vector{Int}`: per-emitter regime
#   in 0..n_regimes (0 = ungroupable, 1 = lowest, n_regimes = highest),
#   in original emitter order.
# - `metadata["mrf_lambda_used"]::Vector{Float64}`: per-group λ actually
#   used (auto or explicit), in `_group_by_dataset` order.
# - `metadata["mrf_regime_means"]::Vector{Vector{Float64}}`: per-group
#   GMM component means in log-density space (sorted ascending), in
#   `_group_by_dataset` order. When `regime_thresholds` is supplied the
#   per-group entry is filled with `NaN`s of length `n_regimes` to signal
#   "manual binning, no GMM fit".
#
# 2D only — DelaunayTriangulation.jl does not provide 3D Voronoi (V7).

import Statistics

"""
    MRFDensityClusterConfig(; n_regimes=2, regime_thresholds=nothing,
                              smoothness_lambda=nothing,
                              graph_kind=:delaunay, graph_k=8,
                              inference=:icm, icm_iters=50,
                              min_points=5, use_3d=false,
                              per_dataset=true, remove_unclustered=false)

Adaptive-density clustering via per-emitter Voronoi density → multi-component
GMM regime assignment → multi-class Potts MRF refinement → connected
components on the foreground.

Designed for SMLM data with multiple density regimes (e.g. one low-density
background + one or more higher-density structures) where a single global
ε would either over-merge or over-split. The MRF smoothness term enforces
spatial coherence — points surrounded by foreground neighbors stay foreground
even if their individual density is borderline (fixes "missing middles");
isolated tight knots in a low-density sea get pulled back to background
(fixes spurious-small-cluster artifacts).

# Fields
- `n_regimes::Int = 2`: number of density regimes. Lowest is treated as
  background/noise; higher regimes form the foreground.
- `regime_thresholds::Union{Nothing, Vector{Float64}} = nothing`: optional
  explicit log-density thresholds, length `n_regimes - 1`, sorted ascending.
  When provided, GMM is bypassed and points are binned directly.
- `smoothness_lambda::Union{Nothing, Float64} = nothing`: MRF smoothness
  weight. When `nothing`, auto-set per group to `max(1e-6, MAD(U_max - U_min))`
  where MAD is the median absolute deviation of the per-emitter unary range.
- `graph_kind::Symbol = :delaunay`: neighbor graph for both the MRF and the
  CC step. `:delaunay` reuses the tessellation from step 1 (free); `:knn`
  builds a symmetrized k-NN graph (uses `graph_k`).
- `graph_k::Int = 8`: k for the kNN graph when `graph_kind=:knn`. Ignored
  for `:delaunay`.
- `inference::Symbol = :icm`: MRF inference algorithm. Only `:icm`
  (Iterated Conditional Modes) is supported in v1; future graph-cut
  inference would land here.
- `icm_iters::Int = 50`: maximum ICM passes; terminates early when no
  point changes label in a full pass.
- `min_points::Int = 5`: minimum cluster size after CC; smaller components
  are demoted to noise.
- `use_3d::Bool = false`: must be `false`. 3D Voronoi tessellation is not
  supported (DelaunayTriangulation.jl is 2D only); passing `true` raises
  `ArgumentError`.
- `per_dataset::Bool = true`: when `true`, run the full pipeline per dataset
  so each cell's density distribution is fit independently.
- `remove_unclustered::Bool = false`: drop noise emitters from the output.

# Outputs
Standard `(smld_out, info)` tuple. Per-cluster cluster ids in `emitter.id`
(`0` = noise / lowest regime / sub-min component, `1..K` = cluster).
HDBSCAN-style metadata stamped onto `smld_out.metadata`:

- `metadata["mrf_regime_per_emitter"]::Vector{Int}` of length `n_locs_in`,
  in original emitter order. Values: 0 (ungroupable / group too small),
  1 (lowest density / background), ..., `n_regimes` (highest density).
- `metadata["mrf_lambda_used"]::Vector{Float64}`: per-group λ used.
- `metadata["mrf_regime_means"]::Vector{Vector{Float64}}`: per-group GMM
  component means (sorted ascending) in log-density space; vector of
  `NaN`s when `regime_thresholds` was provided.

# Example
```julia
# 2-regime auto: GMM finds the foreground/background split per dataset
cfg = MRFDensityClusterConfig()
(smld_out, info) = cluster(smld, cfg)
regimes = smld_out.metadata["mrf_regime_per_emitter"]

# 3-regime with manual thresholds (e.g. learned from training data)
cfg2 = MRFDensityClusterConfig(n_regimes = 3,
                               regime_thresholds = [3.5, 5.0],
                               min_points = 10)
```

See also: [`AbstractClusterConfig`](@ref), [`ClusterInfo`](@ref),
[`cluster`](@ref), [`VoronoiDensityConfig`](@ref) (the read-only sibling
that only exposes per-emitter densities without clustering).
"""
Base.@kwdef struct MRFDensityClusterConfig <: AbstractClusterConfig
    n_regimes::Int = 2
    regime_thresholds::Union{Nothing, Vector{Float64}} = nothing
    smoothness_lambda::Union{Nothing, Float64} = nothing
    graph_kind::Symbol = :delaunay
    graph_k::Int = 8
    inference::Symbol = :icm
    icm_iters::Int = 50
    min_points::Int = 5
    use_3d::Bool = false
    per_dataset::Bool = true
    remove_unclustered::Bool = false
end

# ----------------------------------------------------------------------------
# 1D multi-component Gaussian mixture EM. Returns (means, vars, weights),
# all sorted ascending by mean. Returns `nothing` if the fit fails (variance
# collapses, all points identical density, n < n_components, etc.).
# ----------------------------------------------------------------------------
function _gmm_em_1d(x::Vector{Float64}, k::Int;
                    max_iters::Int = 100, tol::Float64 = 1e-5,
                    var_floor::Float64 = 1e-12)
    n = length(x)
    if n < k || n == 0
        return nothing
    end
    # Bail out if all points identical (variance collapse is unavoidable).
    xmin, xmax = extrema(x)
    if xmin == xmax
        return nothing
    end

    # Init: means from quantiles (deterministic given x, no randomness),
    # variance from overall sample variance, equal weights.
    sx = sort(x)
    μ = Vector{Float64}(undef, k)
    @inbounds for j in 1:k
        q = j / (k + 1)
        idx = clamp(Int(floor(q * n)) + 1, 1, n)
        μ[j] = sx[idx]
    end
    # Overall sample variance.
    mx = sum(x) / n
    σ2_init = sum((xi - mx)^2 for xi in x) / max(n - 1, 1)
    σ2_init = max(σ2_init, var_floor * 100)
    σ2 = fill(σ2_init, k)
    w = fill(1.0 / k, k)

    log_resp = Matrix{Float64}(undef, n, k)
    prev_ll = -Inf
    for iter in 1:max_iters
        # E-step: log responsibilities via logsumexp.
        @inbounds for i in 1:n
            for j in 1:k
                lp = log(max(w[j], 1e-300)) - 0.5 * (log(2π * σ2[j]) + (x[i] - μ[j])^2 / σ2[j])
                log_resp[i, j] = lp
            end
        end
        # logsumexp per row.
        ll = 0.0
        @inbounds for i in 1:n
            row_max = log_resp[i, 1]
            for j in 2:k
                if log_resp[i, j] > row_max
                    row_max = log_resp[i, j]
                end
            end
            s = 0.0
            for j in 1:k
                s += exp(log_resp[i, j] - row_max)
            end
            lse = row_max + log(s)
            ll += lse
            for j in 1:k
                log_resp[i, j] -= lse  # now log of normalized responsibility
            end
        end
        # M-step.
        @inbounds for j in 1:k
            sum_r = 0.0
            sum_rx = 0.0
            for i in 1:n
                r = exp(log_resp[i, j])
                sum_r += r
                sum_rx += r * x[i]
            end
            if sum_r < 1e-12
                # Component starved — bail to caller, who will fall back.
                return nothing
            end
            μ_new = sum_rx / sum_r
            sum_rs = 0.0
            for i in 1:n
                r = exp(log_resp[i, j])
                sum_rs += r * (x[i] - μ_new)^2
            end
            σ2_new = max(sum_rs / sum_r, var_floor)
            μ[j] = μ_new
            σ2[j] = σ2_new
            w[j] = sum_r / n
        end
        # Convergence check.
        if abs(ll - prev_ll) < tol
            break
        end
        prev_ll = ll
    end

    # Sort components ascending by mean so regime 1 = lowest density.
    perm = sortperm(μ)
    μ = μ[perm]
    σ2 = σ2[perm]
    w = w[perm]
    return (means = μ, vars = σ2, weights = w)
end

# ----------------------------------------------------------------------------
# Build the unary cost matrix U[i, k] = -log(w_k * N(x[i] | μ_k, σ_k²)).
# Used in step 3 (MRF refinement) as the data term.
# ----------------------------------------------------------------------------
function _unary_from_gmm(x::Vector{Float64}, fit)
    n = length(x)
    k = length(fit.means)
    U = Matrix{Float64}(undef, n, k)
    @inbounds for j in 1:k
        μj = fit.means[j]
        σ2j = fit.vars[j]
        log_w = log(max(fit.weights[j], 1e-300))
        log_norm = 0.5 * log(2π * σ2j)
        inv_2σ2 = 1.0 / (2.0 * σ2j)
        for i in 1:n
            U[i, j] = log_norm + (x[i] - μj)^2 * inv_2σ2 - log_w
        end
    end
    return U
end

# ----------------------------------------------------------------------------
# Build the unary cost matrix from explicit thresholds. Each emitter goes
# fully into one bin (cost 0) and pays a large penalty for the others.
# This effectively pins step 2's labels; the MRF then only smooths bin
# boundaries via the pairwise term.
# ----------------------------------------------------------------------------
function _unary_from_thresholds(x::Vector{Float64}, thresholds::Vector{Float64}, n_regimes::Int)
    n = length(x)
    U = fill(1.0e6, n, n_regimes)
    @inbounds for i in 1:n
        # Find regime index by binning.
        regime = 1
        for (k, t) in enumerate(thresholds)
            if x[i] >= t
                regime = k + 1
            end
        end
        U[i, regime] = 0.0
    end
    return U
end

# ----------------------------------------------------------------------------
# Build a Vector{Vector{Int}} adjacency list for n nodes given an iterable
# of (i, j) edges with i < j.
# ----------------------------------------------------------------------------
function _build_adj(n::Int, edges)
    adj = [Int[] for _ in 1:n]
    @inbounds for (i, j) in edges
        push!(adj[i], j)
        push!(adj[j], i)
    end
    return adj
end

# ----------------------------------------------------------------------------
# Delaunay edges for the given triangulation, deduplicated to (i<j).
# ----------------------------------------------------------------------------
function _delaunay_edges(tri)
    edges = Tuple{Int, Int}[]
    for e in DelaunayTriangulation.each_solid_edge(tri)
        i = min(e[1], e[2])
        j = max(e[1], e[2])
        i == j && continue
        push!(edges, (i, j))
    end
    # Dedup (each_solid_edge may yield each edge once or twice depending on
    # version; canonicalize defensively).
    sort!(edges)
    unique!(edges)
    return edges
end

# ----------------------------------------------------------------------------
# Symmetrized k-NN edge list from a 2×n coordinate matrix. Excludes self.
# ----------------------------------------------------------------------------
function _knn_edges(X::Matrix{Float64}, k::Int)
    n = size(X, 2)
    k_use = min(k, n - 1)
    if k_use < 1
        return Tuple{Int, Int}[]
    end
    tree = NearestNeighbors.KDTree(X)
    edges = Set{Tuple{Int, Int}}()
    @inbounds for i in 1:n
        idxs, _ = NearestNeighbors.knn(tree, view(X, :, i), k_use + 1, true)
        for jj in idxs
            jj == i && continue
            a = min(i, jj)
            b = max(i, jj)
            push!(edges, (a, b))
        end
    end
    return collect(edges)
end

# ----------------------------------------------------------------------------
# Iterated Conditional Modes for multi-class Potts MRF. Initialize labels
# from argmin of unary; iterate, updating each point's label to argmin of
# (U[i, k] + λ * #neighbors_with_label_not_k). Terminate when no changes
# in a full pass or after `max_iters`.
# ----------------------------------------------------------------------------
function _icm_potts!(labels::Vector{Int}, U::Matrix{Float64},
                     adj::Vector{Vector{Int}}, λ::Float64, max_iters::Int)
    n, k = size(U)
    @inbounds for i in 1:n
        # Initialize from argmin unary.
        best_j = 1
        best_v = U[i, 1]
        for j in 2:k
            if U[i, j] < best_v
                best_v = U[i, j]
                best_j = j
            end
        end
        labels[i] = best_j
    end
    for _ in 1:max_iters
        changed = 0
        @inbounds for i in 1:n
            cur = labels[i]
            best_j = cur
            # Cost of staying = U[i, cur] + λ * count of neighbors with label != cur.
            best_v = Inf
            for j in 1:k
                neigh_diff = 0
                for q in adj[i]
                    if labels[q] != j
                        neigh_diff += 1
                    end
                end
                v = U[i, j] + λ * neigh_diff
                if v < best_v
                    best_v = v
                    best_j = j
                end
            end
            if best_j != cur
                labels[i] = best_j
                changed += 1
            end
        end
        changed == 0 && break
    end
    return labels
end

# ----------------------------------------------------------------------------
# Connected components via BFS on `adj`, restricted to nodes where
# `is_foreground[i]` is true (and edges to other foreground nodes).
# Returns (raw_labels, n_components) where raw_labels[i] = 0 if background
# or 1..n_components for foreground.
# ----------------------------------------------------------------------------
function _connected_components_bfs(adj::Vector{Vector{Int}}, is_foreground::AbstractVector{Bool})
    n = length(adj)
    raw_labels = zeros(Int, n)
    n_comp = 0
    queue = Int[]
    @inbounds for seed in 1:n
        (is_foreground[seed] && raw_labels[seed] == 0) || continue
        n_comp += 1
        raw_labels[seed] = n_comp
        empty!(queue)
        push!(queue, seed)
        while !isempty(queue)
            v = pop!(queue)
            for w in adj[v]
                if is_foreground[w] && raw_labels[w] == 0
                    raw_labels[w] = n_comp
                    push!(queue, w)
                end
            end
        end
    end
    return raw_labels, n_comp
end

# ----------------------------------------------------------------------------
# SMLD-facing dispatch.
# ----------------------------------------------------------------------------
function cluster(smld::SMLMData.BasicSMLD, cfg::MRFDensityClusterConfig)
    t0 = time_ns()
    smld = SMLMData.BasicSMLD(deepcopy(smld.emitters), smld.camera,
                              smld.n_frames, smld.n_datasets,
                              deepcopy(smld.metadata))
    n_in = length(smld.emitters)

    # Argument validation.
    cfg.n_regimes >= 2 ||
        throw(ArgumentError("MRFDensityClusterConfig.n_regimes must be ≥ 2 (got $(cfg.n_regimes))"))
    cfg.graph_kind in (:delaunay, :knn) ||
        throw(ArgumentError("MRFDensityClusterConfig.graph_kind must be :delaunay or :knn (got $(cfg.graph_kind))"))
    cfg.inference === :icm ||
        throw(ArgumentError("MRFDensityClusterConfig.inference: only :icm currently supported (got $(cfg.inference))"))
    cfg.graph_k >= 1 ||
        throw(ArgumentError("MRFDensityClusterConfig.graph_k must be ≥ 1 (got $(cfg.graph_k))"))
    cfg.icm_iters >= 1 ||
        throw(ArgumentError("MRFDensityClusterConfig.icm_iters must be ≥ 1 (got $(cfg.icm_iters))"))
    cfg.min_points >= 1 ||
        throw(ArgumentError("MRFDensityClusterConfig.min_points must be ≥ 1 (got $(cfg.min_points))"))
    cfg.use_3d &&
        throw(ArgumentError(
            "MRFDensityClusterConfig does not support use_3d=true. " *
            "Voronoi tessellation is 2D only (V7); use DBSCANConfig or " *
            "HierarchicalConfig with use_3d=true for 3D data."))
    if cfg.smoothness_lambda !== nothing && !(cfg.smoothness_lambda > 0)
        throw(ArgumentError("MRFDensityClusterConfig.smoothness_lambda must be > 0 when set (got $(cfg.smoothness_lambda))"))
    end
    if cfg.regime_thresholds !== nothing
        ts = cfg.regime_thresholds
        length(ts) == cfg.n_regimes - 1 ||
            throw(ArgumentError("MRFDensityClusterConfig.regime_thresholds length must be n_regimes - 1 = $(cfg.n_regimes - 1) (got $(length(ts)))"))
        issorted(ts) ||
            throw(ArgumentError("MRFDensityClusterConfig.regime_thresholds must be sorted ascending (got $ts)"))
    end

    # Per-emitter regime, in original emitter order. 0 = ungroupable.
    regime_per_emitter = zeros(Int, n_in)
    lambda_used_per_group = Float64[]
    means_per_group = Vector{Vector{Float64}}()

    groups = _group_by_dataset(smld, cfg.per_dataset)

    cluster_sizes = Int[]
    n_clustered = 0

    for idxs in groups
        n = length(idxs)
        if n < 3
            # Too small for tessellation — all noise (regime 0 is the default
            # already from `zeros(Int, n_in)` above).
            @inbounds for i in idxs
                smld.emitters[i].id = 0
            end
            push!(lambda_used_per_group, NaN)
            push!(means_per_group, fill(NaN, cfg.n_regimes))
            continue
        end

        sub = view(smld.emitters, idxs)
        # Step 1: per-emitter Voronoi areas (and the triangulation, kept for
        # the Delaunay edge graph). Raises ArgumentError on duplicate (x,y).
        areas, tri = _voronoi_areas(sub)

        # log densities (μm⁻²) for all-finite-positive areas; NaN for invalid.
        log_rho = Vector{Float64}(undef, n)
        valid_mask = falses(n)
        @inbounds for j in 1:n
            a = areas[j]
            if isfinite(a) && a > 0
                log_rho[j] = log(1.0 / a)
                valid_mask[j] = true
            else
                log_rho[j] = NaN
            end
        end
        n_valid = count(valid_mask)
        if n_valid < cfg.n_regimes
            # Not enough valid points to fit the regime structure — treat all as noise.
            @inbounds for i in idxs
                smld.emitters[i].id = 0
            end
            push!(lambda_used_per_group, NaN)
            push!(means_per_group, fill(NaN, cfg.n_regimes))
            continue
        end

        # Step 2: regime assignment. Build U over the valid emitters.
        valid_idxs = findall(valid_mask)  # local indices 1..n
        log_rho_valid = log_rho[valid_idxs]
        local fit_means_for_metadata::Vector{Float64}
        local U::Matrix{Float64}
        if cfg.regime_thresholds !== nothing
            U = _unary_from_thresholds(log_rho_valid, cfg.regime_thresholds, cfg.n_regimes)
            fit_means_for_metadata = fill(NaN, cfg.n_regimes)
        else
            fit = _gmm_em_1d(log_rho_valid, cfg.n_regimes)
            if fit === nothing
                # GMM failed → fall back: all valid points to regime 1 (lowest).
                @inbounds for j in valid_idxs
                    regime_per_emitter[idxs[j]] = 1
                end
                @inbounds for i in idxs
                    smld.emitters[i].id = 0
                end
                push!(lambda_used_per_group, NaN)
                push!(means_per_group, fill(NaN, cfg.n_regimes))
                continue
            end
            U = _unary_from_gmm(log_rho_valid, fit)
            fit_means_for_metadata = collect(fit.means)
        end

        # Step 3: build neighbor graph (over the FULL n points; the MRF/CC
        # only operates on valid_idxs in practice but the adjacency list
        # carries everyone for index simplicity).
        edges = if cfg.graph_kind === :delaunay
            _delaunay_edges(tri)
        else
            X_local = _coords_matrix(sub, false)
            _knn_edges(X_local, cfg.graph_k)
        end
        adj = _build_adj(n, edges)

        # Filter the adjacency to only valid nodes for the MRF — but ICM is
        # cleaner if we run it on a re-indexed valid-only subgraph. Build
        # local mapping (full index → valid index, or 0).
        full_to_valid = zeros(Int, n)
        @inbounds for (vidx, j) in enumerate(valid_idxs)
            full_to_valid[j] = vidx
        end
        n_v = length(valid_idxs)
        adj_valid = [Int[] for _ in 1:n_v]
        @inbounds for j in 1:n
            full_to_valid[j] == 0 && continue
            for w in adj[j]
                full_to_valid[w] == 0 && continue
                push!(adj_valid[full_to_valid[j]], full_to_valid[w])
            end
        end

        # Auto-λ heuristic: MAD of per-emitter unary range over valid points.
        λ = if cfg.smoothness_lambda === nothing
            ranges = Vector{Float64}(undef, n_v)
            @inbounds for i in 1:n_v
                vmin = U[i, 1]; vmax = U[i, 1]
                for k in 2:cfg.n_regimes
                    v = U[i, k]
                    if v < vmin; vmin = v; end
                    if v > vmax; vmax = v; end
                end
                ranges[i] = vmax - vmin
            end
            med = Statistics.median(ranges)
            mad = Statistics.median(abs.(ranges .- med))
            max(1e-6, mad)
        else
            cfg.smoothness_lambda
        end
        push!(lambda_used_per_group, λ)
        push!(means_per_group, fit_means_for_metadata)

        # Step 3 (cont.): ICM.
        labels_valid = Vector{Int}(undef, n_v)
        _icm_potts!(labels_valid, U, adj_valid, λ, cfg.icm_iters)

        # Stamp regimes back to per-emitter (original emitter order).
        @inbounds for (vidx, j) in enumerate(valid_idxs)
            regime_per_emitter[idxs[j]] = labels_valid[vidx]
        end

        # Step 4: connected components on foreground (regime ≥ 2) over the
        # full-group neighbor graph (so foreground points connected via a
        # short path through other foreground points get merged correctly).
        is_foreground = falses(n)
        @inbounds for (vidx, j) in enumerate(valid_idxs)
            if labels_valid[vidx] >= 2
                is_foreground[j] = true
            end
        end
        raw_labels, n_raw = _connected_components_bfs(adj, is_foreground)

        # Apply min_points size filter and compact-relabel.
        raw_counts = zeros(Int, n_raw)
        @inbounds for l in raw_labels
            l > 0 && (raw_counts[l] += 1)
        end
        label_map, added = _compact_relabel!(cluster_sizes, raw_counts, cfg.min_points)
        n_clustered += added

        @inbounds for (j, i) in pairs(idxs)
            rl = raw_labels[j]
            smld.emitters[i].id = rl > 0 ? label_map[rl] : 0
        end
    end

    n_clusters = length(cluster_sizes)
    n_noise = n_in - n_clustered
    smld_out = _build_output(smld, cfg.remove_unclustered)
    smld_out.metadata["mrf_regime_per_emitter"] = regime_per_emitter
    smld_out.metadata["mrf_lambda_used"]        = lambda_used_per_group
    smld_out.metadata["mrf_regime_means"]       = means_per_group

    info = ClusterInfo(
        n_in,
        n_clustered,
        n_noise,
        n_clusters,
        cluster_sizes,
        :mrf_density,
        (time_ns() - t0) / 1e9,
    )
    return smld_out, info
end
