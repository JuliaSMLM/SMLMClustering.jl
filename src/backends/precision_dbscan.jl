# Precision-weighted (σ-aware) DBSCAN backend.
#
# Two layers:
#
#   1. A reusable, geometry-only neighbor cache — `PrecisionNeighborGraph` +
#      `build_precision_neighbor_graph` — and an exact re-threshold+label pass
#      (`precision_dbscan_labels`). The graph caches candidate pairs (i, j, raw
#      Euclidean distance d) once, up to a caller-supplied `max_radius`; the label
#      pass re-thresholds `d < nsigma·(σ_eff_i + σ_eff_j)` on the cache and never
#      rebuilds the tree. This is the primitive SMLMBaGoL reuses across the ~10
#      E-steps of its τ-finder (only σ_eff/nsigma vary; coordinates are fixed).
#
#   2. An idiomatic `cluster(smld, ::PrecisionDBSCANConfig)` wrapper for the
#      SMLD-facing path: derives σ_eff from each emitter's localization precision,
#      builds the graph internally, and writes cluster ids onto `emitter.id`.
#
# Unlike the Euclidean `DBSCANConfig` (which wraps `Clustering.dbscan`), the metric
# here is precision-weighted, so this is a self-contained implementation. The
# neighbor prepass is threaded (mirrors `local_contrast.jl`); the label pass is
# deterministic: union-find connected components for `min_points == 0` (order-free,
# bit-exact regardless of thread scheduling) and a deterministic core-point
# variant for `min_points ≥ 1`.

import Base.Threads: @threads

"""
    PrecisionNeighborGraph

Geometry-only neighbor cache for precision-weighted DBSCAN, built once by
[`build_precision_neighbor_graph`](@ref) and reused across many
[`precision_dbscan_labels`](@ref) calls with varying `σ_eff` / `nsigma`.

Treat it as an **opaque handle**: build it with `build_precision_neighbor_graph`
and pass it back to `precision_dbscan_labels`. Its fields (`n`, `offsets`,
`neighbors`, `dists`, `max_radius`, `dims`) are an **implementation detail** — a CSR
store of every candidate pair with raw Euclidean `d ≤ max_radius` (full adjacency,
both directions) — and **may change between releases**; you should not need to read them.

The cache is a valid superset for a label pass iff
`nsigma·(σ_eff_i + σ_eff_j) ≤ max_radius` for every pair; the sufficient rule is
`max_radius ≥ nsigma · 2 · maximum(σ_eff)` (checked at label time when
`check_superset = true`).
"""
struct PrecisionNeighborGraph
    n::Int
    offsets::Vector{Int}
    neighbors::Vector{Int}
    dists::Vector{Float64}
    max_radius::Float64
    dims::Int
end

function Base.show(io::IO, g::PrecisionNeighborGraph)
    nedges = length(g.neighbors) ÷ 2
    print(io, "PrecisionNeighborGraph(", g.n, " pts, ", nedges,
          " undirected pairs ≤ ", round(g.max_radius, sigdigits = 4),
          ", ", g.dims, "D)")
end

"""
    build_precision_neighbor_graph(coords::AbstractMatrix{<:Real}, max_radius::Real)
        -> PrecisionNeighborGraph

Build the geometry-only neighbor cache. `coords` is a `D×N` matrix (columns are
points, `D ∈ {2, 3}`) in whatever length unit you will also use for `max_radius`
and `σ_eff`. Every pair within `max_radius` (raw Euclidean) is cached with its
distance; the σ-weighted threshold is applied later, per call, by
[`precision_dbscan_labels`](@ref).

The per-point range queries are the expensive part and are run in parallel
(`@threads`); each thread writes only its own point's neighbor list, so no
per-thread scratch is needed and the result is independent of thread scheduling.
Neighbor lists are sorted by index, making the CSR — and hence the labels —
reproducible.
"""
function build_precision_neighbor_graph(coords::AbstractMatrix{<:Real}, max_radius::Real)
    D, n = size(coords)
    (D == 2 || D == 3) ||
        throw(ArgumentError("coords must be 2×N or 3×N (columns = points); got $(D)×$(n)"))
    max_radius > 0 || throw(ArgumentError("max_radius must be > 0; got $max_radius"))
    r = Float64(max_radius)
    if n == 0
        return PrecisionNeighborGraph(0, Int[1], Int[], Float64[], r, D)
    end
    X = Matrix{Float64}(coords)          # dense Float64 D×N for the KDTree
    tree = NearestNeighbors.KDTree(X)
    nbr_lists = Vector{Vector{Int}}(undef, n)
    dst_lists = Vector{Vector{Float64}}(undef, n)
    @threads for i in 1:n
        idx = NearestNeighbors.inrange(tree, view(X, :, i), r)   # includes self
        js = Int[]
        for j in idx
            j == i || push!(js, j)
        end
        sort!(js)                        # deterministic CSR order
        ds = Vector{Float64}(undef, length(js))
        @inbounds for (t, j) in enumerate(js)
            acc = 0.0
            for d in 1:D
                δ = X[d, i] - X[d, j]
                acc += δ * δ
            end
            ds[t] = sqrt(acc)
        end
        nbr_lists[i] = js
        dst_lists[i] = ds
    end
    # Flatten per-point lists into CSR (serial, O(total edges) — cheap next to the
    # threaded range queries above).
    offsets = Vector{Int}(undef, n + 1)
    offsets[1] = 1
    @inbounds for i in 1:n
        offsets[i + 1] = offsets[i] + length(nbr_lists[i])
    end
    total = offsets[n + 1] - 1
    neighbors = Vector{Int}(undef, total)
    dists = Vector{Float64}(undef, total)
    @inbounds for i in 1:n
        base = offsets[i] - 1
        js = nbr_lists[i]; ds = dst_lists[i]
        for t in eachindex(js)
            neighbors[base + t] = js[t]
            dists[base + t] = ds[t]
        end
    end
    return PrecisionNeighborGraph(n, offsets, neighbors, dists, r, D)
end

# ---- union-find (path halving + union by size) ------------------------------

@inline function _uf_find(parent::Vector{Int}, x::Int)
    @inbounds while parent[x] != x
        parent[x] = parent[parent[x]]
        x = parent[x]
    end
    return x
end

@inline function _uf_union!(parent::Vector{Int}, sz::Vector{Int}, a::Int, b::Int)
    ra = _uf_find(parent, a); rb = _uf_find(parent, b)
    ra == rb && return nothing
    @inbounds if sz[ra] < sz[rb]
        ra, rb = rb, ra
    end
    @inbounds parent[rb] = ra
    @inbounds sz[ra] += sz[rb]
    return nothing
end

"""
    precision_dbscan_labels(g::PrecisionNeighborGraph, σ_eff, nsigma;
                            min_points = 1, check_superset = true) -> Vector{Int}
    precision_dbscan_labels!(labels, g, σ_eff, nsigma; min_points = 1, check_superset = true)

Label the points of a prebuilt [`PrecisionNeighborGraph`](@ref) by re-thresholding
the cached pairs: an edge `(i, j)` is **active** iff
`g.dists < nsigma·(σ_eff[i] + σ_eff[j])`. `σ_eff` (length `g.n`) and `nsigma` may
change on every call **without rebuilding `g`** — that reuse is the point of the
primitive.

- `min_points == 0`: connected components of the active graph, via union-find. The
  partition is **order-free / bit-exact** — independent of thread scheduling and
  edge order. Every point is labeled `1..K` (a singleton is its own component; no
  noise).
- `min_points ≥ 1`: core-point DBSCAN with the classical (self-inclusive) `minPts`,
  identical to [`DBSCANConfig`](@ref). A point is *core* iff its neighborhood — itself
  plus its active neighbors — has `≥ min_points` members (i.e. active degree
  `≥ min_points − 1`); core points sharing an active edge merge; a non-core *border*
  point joins the **lowest-id** adjacent core cluster (deterministic tie-break); a
  non-core point with no active core neighbor is noise (`0`).

Cluster ids are canonical: `1..K` in ascending order of first appearance by point
index, so the integer labels themselves are reproducible.

`check_superset` (default `true`) asserts `nsigma·2·maximum(σ_eff) ≤ g.max_radius`
and errors if the cache is too tight to contain every active pair (see
[`PrecisionNeighborGraph`](@ref)). The mutating form fills a caller-supplied
`labels::Vector{Int}` of length `g.n`.
"""
function precision_dbscan_labels(g::PrecisionNeighborGraph, σ_eff::AbstractVector{<:Real},
                                 nsigma::Real; min_points::Int = 1, check_superset::Bool = true)
    labels = Vector{Int}(undef, g.n)
    return precision_dbscan_labels!(labels, g, σ_eff, nsigma;
                                    min_points = min_points, check_superset = check_superset)
end

"""
    precision_dbscan_labels!(labels, g, σ_eff, nsigma; min_points=1, check_superset=true) -> labels

In-place form of [`precision_dbscan_labels`](@ref): fills the caller-supplied
`labels::Vector{Int}` (which must have length `g.n`) and returns it. Use this in a
hot loop that reuses one [`PrecisionNeighborGraph`](@ref) across many calls to avoid
reallocating the label vector each time.
"""
function precision_dbscan_labels!(labels::Vector{Int}, g::PrecisionNeighborGraph,
                                  σ_eff::AbstractVector{<:Real}, nsigma::Real;
                                  min_points::Int = 1, check_superset::Bool = true)
    n = g.n
    length(labels) == n ||
        throw(ArgumentError("labels length $(length(labels)) ≠ graph n $n"))
    length(σ_eff) == n ||
        throw(ArgumentError("σ_eff length $(length(σ_eff)) ≠ graph n $n"))
    nsigma > 0 || throw(ArgumentError("nsigma must be > 0; got $nsigma"))
    min_points >= 0 || throw(ArgumentError("min_points must be ≥ 0; got $min_points"))
    n == 0 && return labels
    if check_superset
        # The cache holds every pair with raw d ≤ max_radius; a pair is active iff
        # d < nsigma·(σ_i+σ_j) ≤ nsigma·2·maximum(σ_eff). So the cache is a valid
        # superset iff that bound ≤ max_radius. The 1e-9 is float slack for callers
        # (e.g. the config wrapper) that build max_radius == need with a different op order.
        need = Float64(nsigma) * 2 * maximum(σ_eff)
        need <= g.max_radius * (1 + 1e-9) || throw(ArgumentError(
            "precision_dbscan_labels: cache too tight — nsigma·2·maximum(σ_eff) = " *
            "$need exceeds graph max_radius = $(g.max_radius); rebuild the graph with " *
            "max_radius ≥ $need (a pair could be active yet uncached)."))
    end
    ns = Float64(nsigma)
    off = g.offsets; nb = g.neighbors; ds = g.dists
    parent = collect(1:n)
    sz = ones(Int, n)

    if min_points <= 0
        # Connected components over active edges (each undirected edge once, j > i).
        @inbounds for i in 1:n
            σi = Float64(σ_eff[i])
            for t in off[i]:(off[i + 1] - 1)
                j = nb[t]
                j > i || continue
                ds[t] < ns * (σi + Float64(σ_eff[j])) && _uf_union!(parent, sz, i, j)
            end
        end
        return _canonicalize_all!(labels, parent, n)
    end

    # Core points: classical DBSCAN minPts — a point is core iff its neighborhood
    # (itself + active neighbors) has ≥ min_points members, i.e. active degree ≥
    # min_points - 1. Matches DBSCANConfig / Clustering.dbscan (self-inclusive count).
    core = falses(n)
    @inbounds for i in 1:n
        σi = Float64(σ_eff[i]); c = 0
        for t in off[i]:(off[i + 1] - 1)
            ds[t] < ns * (σi + Float64(σ_eff[nb[t]])) && (c += 1)
        end
        core[i] = c >= min_points - 1
    end
    # Merge core points connected by an active edge.
    @inbounds for i in 1:n
        core[i] || continue
        σi = Float64(σ_eff[i])
        for t in off[i]:(off[i + 1] - 1)
            j = nb[t]
            j > i || continue
            (core[j] && ds[t] < ns * (σi + Float64(σ_eff[j]))) && _uf_union!(parent, sz, i, j)
        end
    end
    # Canonical labels for core components (ascending first-appearance).
    fill!(labels, 0)
    rootlabel = zeros(Int, n); nextlab = 0
    @inbounds for i in 1:n
        core[i] || continue
        r = _uf_find(parent, i)
        rootlabel[r] == 0 && (nextlab += 1; rootlabel[r] = nextlab)
        labels[i] = rootlabel[r]
    end
    # Border points → lowest-id adjacent active core cluster; else noise (0).
    @inbounds for i in 1:n
        core[i] && continue
        σi = Float64(σ_eff[i]); best = 0
        for t in off[i]:(off[i + 1] - 1)
            j = nb[t]
            (core[j] && ds[t] < ns * (σi + Float64(σ_eff[j]))) || continue
            lj = labels[j]
            (best == 0 || lj < best) && (best = lj)
        end
        labels[i] = best
    end
    return labels
end

# Relabel union-find roots to 1..K by ascending first-appearance (min_points==0).
function _canonicalize_all!(labels::Vector{Int}, parent::Vector{Int}, n::Int)
    rootlabel = zeros(Int, n); nextlab = 0
    @inbounds for i in 1:n
        r = _uf_find(parent, i)
        rootlabel[r] == 0 && (nextlab += 1; rootlabel[r] = nextlab)
        labels[i] = rootlabel[r]
    end
    return labels
end

# ---- SMLD-facing config + cluster() -----------------------------------------

"""
    PrecisionDBSCANConfig(; nsigma, min_points=5, use_3d=false, per_dataset=true,
                          remove_unclustered=false)

Configuration for **precision-weighted DBSCAN** of SMLM localizations: the
neighbor test is `‖pᵢ − pⱼ‖ < nsigma · (σ_effᵢ + σ_effⱼ)`, where each emitter's
`σ_eff` is the geometric mean of its per-axis localization precisions
(`√(σ_x·σ_y)` in 2D, `∛(σ_x·σ_y·σ_z)` in 3D). Unlike [`DBSCANConfig`](@ref)'s
fixed `eps_nm`, the neighborhood adapts to each localization's uncertainty.

# Fields
- `nsigma::Float64`: neighbor radius in units of the summed precision `σ_effᵢ + σ_effⱼ`.
- `min_points::Int = 5`: classical DBSCAN `minPts` — the minimum neighborhood size
  (the point itself plus its active neighbors) for a core point — and the minimum
  cluster size (clusters smaller than `min_points` become noise). Identical in meaning
  and default to [`DBSCANConfig`](@ref).
- `use_3d::Bool = false`: cluster in (x, y, z) using `σ_z`; requires `Emitter3DFit`.
- `per_dataset::Bool = true`: cluster within each `dataset` index independently.
- `remove_unclustered::Bool = false`: drop noise emitters (`id == 0`) from the output.

Cluster ids are written to `emitter.id` (`0` = noise, `1..K` = clusters); returns
`(smld_out, ClusterInfo)` with `algorithm = :precision_dbscan`.

For the lower-level reuse-the-graph primitive (build once, relabel many times with
varying `σ_eff`/`nsigma`), use [`build_precision_neighbor_graph`](@ref) +
[`precision_dbscan_labels`](@ref) directly.

See also: [`DBSCANConfig`](@ref), [`AbstractClusterConfig`](@ref), [`cluster`](@ref).
"""
Base.@kwdef struct PrecisionDBSCANConfig <: AbstractClusterConfig
    nsigma::Float64
    min_points::Int = 5
    use_3d::Bool = false
    per_dataset::Bool = true
    remove_unclustered::Bool = false
end

# Per-emitter σ_eff: geometric mean of per-axis localization precisions (μm).
function _sigma_eff(emitters::AbstractVector{<:SMLMData.AbstractEmitter}, use_3d::Bool)
    n = length(emitters)
    s = Vector{Float64}(undef, n)
    if use_3d
        isempty(emitters) || hasproperty(first(emitters), :σ_z) ||
            error("use_3d=true requires 3D emitters with σ_z (e.g. Emitter3DFit); got $(eltype(emitters)).")
        @inbounds for i in 1:n
            e = emitters[i]
            s[i] = cbrt(e.σ_x * e.σ_y * e.σ_z)
        end
    else
        @inbounds for i in 1:n
            e = emitters[i]
            s[i] = sqrt(e.σ_x * e.σ_y)
        end
    end
    return s
end

function cluster(smld::SMLMData.BasicSMLD, cfg::PrecisionDBSCANConfig)
    t0 = time_ns()
    # Non-mutating: labels go onto a fresh copy, not the caller's input.
    smld = SMLMData.BasicSMLD(deepcopy(smld.emitters), smld.camera,
                              smld.n_frames, smld.n_datasets, smld.metadata)
    n_in = length(smld.emitters)
    cfg.nsigma > 0 || throw(ArgumentError("PrecisionDBSCANConfig.nsigma must be > 0 (got $(cfg.nsigma))"))
    cfg.min_points >= 1 ||
        throw(ArgumentError("PrecisionDBSCANConfig.min_points must be ≥ 1 (got $(cfg.min_points))"))

    groups = _group_by_dataset(smld, cfg.per_dataset)
    cluster_sizes = Int[]
    n_clustered = 0

    for idxs in groups
        isempty(idxs) && continue
        sub = view(smld.emitters, idxs)
        coords = _coords_matrix(sub, cfg.use_3d)
        σ_eff = _sigma_eff(sub, cfg.use_3d)
        # Build the cache exactly large enough for this call's threshold, then label.
        max_radius = cfg.nsigma * 2 * maximum(σ_eff)
        max_radius > 0 ||
            throw(ArgumentError("PrecisionDBSCANConfig: all localization precisions are zero; " *
                                "cannot form a precision-weighted neighborhood."))
        g = build_precision_neighbor_graph(coords, max_radius)
        labels = precision_dbscan_labels(g, σ_eff, cfg.nsigma;
                                         min_points = cfg.min_points, check_superset = true)

        # Drop clusters below min_points (mirror DBSCANConfig), compact-relabel,
        # and write ids onto emitter.id. Sizes recounted from the labels.
        k_raw = maximum(labels; init = 0)
        group_sizes = zeros(Int, k_raw)
        @inbounds for a in labels
            a == 0 && continue
            group_sizes[a] += 1
        end
        label_map, added = _compact_relabel!(cluster_sizes, group_sizes, cfg.min_points)
        @inbounds for (j, i) in pairs(idxs)
            a = labels[j]
            smld.emitters[i].id = a == 0 ? 0 : label_map[a]
        end
        n_clustered += added
    end

    n_clusters = length(cluster_sizes)
    n_noise = n_in - n_clustered
    smld_out = _build_output(smld, cfg.remove_unclustered)
    info = ClusterInfo(n_in, n_clustered, n_noise, n_clusters, cluster_sizes,
                       :precision_dbscan, (time_ns() - t0) / 1e9)
    return smld_out, info
end
