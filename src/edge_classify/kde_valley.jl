"""
Internal: `METHOD = "kde_valley"` — the validated adaptive edge gate promoted
from the paper-genmab-hexabody production wrapper (`_analyze_edge_kde`).

Three stages, in the validated order (inverts the v1 reflect-then-gate order —
the KDE gate runs on the ORIGINAL cloud, pre-reflection):

1. **KDE-valley density gate** — continuous Gaussian-KDE density on the original
   cloud, thresholded at the background/cell valley (left-base of the cell mode
   in the log-density histogram). Per-FOV adaptive: handles the ~6× MAP-N density
   spread across FOVs with no per-cell tuning, where a single global `RHO_K_THRESH`
   erases sparse cells while passing dense ones.
2. **Footprint fill** — rasterize the kept tissue, morphological-close thin necks,
   flood-fill enclosed holes. Seals low-density channels so the alpha-shape does
   not leak into interior voids.
3. **Enclosure reclass** — a background point enclosed by the cell (≥
   `ENCLOSURE_MIN_HITS` of 8 rays hit cell tissue before the field edge) is folded
   into `interior`.

The boundary/membrane geometry (FOV reflection, alpha-shape, point-in-polygon,
membrane band) is reused unchanged: stages 1–2 select a footprint SUBSET, the v1
outer-polygon classifier runs on that subset (internal k-NN gate disabled via
`RHO_K_THRESH = 0`, matching the wrapper), and stage 3 post-processes the full set.

Topology contract (i): `class` is authoritative and folds in the enclosure-recovered
interiors; `inside_outer`/`dist_to_outer_um` stay strictly geometric (inside the
alpha outer loop). The enclosure-recovered set is exactly
`class == "interior" && inside_outer == false`. `in_cell` carries topological
cell membership (geometric-inside ∪ enclosure-recovered).

Defaults reproduce the A431 dSTORM validated set (commit 45b0690). dSTORM path
only; DNA-PAINT uses `outer_polygon` with a per-FOV density quantile, not this.
"""

# ---- KDE density (continuous Gaussian, per-point query) -----------------------

# Continuous Gaussian-KDE density at each point (fixed bandwidth `sigma`, µm).
# Per-point range query so it is memory-safe on dense clouds (~290k pts) and
# unbiased at edges/clusters (unlike k-NN ρ_K). Self-contribution subtracted.
function _kde_density(X::AbstractMatrix{Float64}, tree, sigma::Float64;
                      rmax_sigma::Float64 = 3.0)
    n = size(X, 2)
    rho = zeros(n)
    inv2s2 = 1 / (2 * sigma^2)
    rmax = rmax_sigma * sigma
    pt = zeros(2)
    for i in 1:n
        pt[1] = X[1, i]; pt[2] = X[2, i]
        idx = NearestNeighbors.inrange(tree, pt, rmax)
        s = 0.0
        @inbounds for j in idx
            d2 = (X[1, i] - X[1, j])^2 + (X[2, i] - X[2, j])^2
            s += exp(-d2 * inv2s2)
        end
        rho[i] = (s - 1.0) / (2 * pi * sigma^2)   # subtract self
    end
    return rho
end

# ---- Valley threshold ---------------------------------------------------------

# Left base of the dominant (cell) mode: scan left from the global peak of the
# smoothed log-density histogram to the first bin below `floorfrac` of the peak
# (the background/cell gap). Robust where background is a small class.
function _kde_leftbase(vals::AbstractVector{<:Real}; nbins::Int = 140,
                       floorfrac::Float64 = 0.05, smooth::Int = 4)
    lo, hi = extrema(vals)
    lo == hi && return lo
    edges = collect(range(lo, hi, length = nbins + 1))
    centers = (edges[1:end-1] .+ edges[2:end]) ./ 2
    h = zeros(Int, nbins)
    for v in vals
        b = clamp(searchsortedlast(edges, v), 1, nbins)
        h[b] += 1
    end
    hs = [Statistics.mean(@view h[max(1, i - smooth):min(nbins, i + smooth)]) for i in 1:nbins]
    pk = argmax(hs)
    thr = floorfrac * hs[pk]
    i = pk
    while i > 1 && hs[i] >= thr
        i -= 1
    end
    return centers[i]
end

# ---- Binary morphology + footprint fill (occupancy grid, no deps) -------------

# Binary dilate (`dilate=true`) / erode (`dilate=false`) on an occupancy grid by
# a square structuring element of radius `r`.
function _morph(g::AbstractMatrix{Bool}, r::Int, dilate::Bool)
    nx, ny = size(g)
    out = falses(nx, ny)
    for i in 1:nx, j in 1:ny
        g[i, j] || continue
        if dilate
            for di in -r:r, dj in -r:r
                ii = i + di; jj = j + dj
                (1 <= ii <= nx && 1 <= jj <= ny) && (out[ii, jj] = true)
            end
        else
            ok = true
            for di in -r:r, dj in -r:r
                ii = i + di; jj = j + dj
                (!(1 <= ii <= nx && 1 <= jj <= ny) || !g[ii, jj]) && (ok = false)
            end
            out[i, j] = ok
        end
    end
    return out
end

# Footprint of the kept tissue: rasterize → dilate-seal thin necks → flood-fill
# enclosed holes from the grid border (original boundary preserved, no erosion).
# Returns a per-point in-footprint mask aligned to (xs, ys).
function _footprint_fill(xs::Vector{Float64}, ys::Vector{Float64},
                         tissue::AbstractVector{Bool}; bin::Float64 = 0.2,
                         closing::Int = 3)
    x0, x1 = extrema(xs); y0, y1 = extrema(ys)
    nx = max(1, ceil(Int, (x1 - x0) / bin) + 1)
    ny = max(1, ceil(Int, (y1 - y0) / bin) + 1)
    bx(x) = clamp(floor(Int, (x - x0) / bin) + 1, 1, nx)
    by(y) = clamp(floor(Int, (y - y0) / bin) + 1, 1, ny)
    occ = falses(nx, ny)
    @inbounds for i in eachindex(xs)
        tissue[i] && (occ[bx(xs[i]), by(ys[i])] = true)
    end
    dil = _morph(occ, closing, true)
    seen = falses(nx, ny)
    stk = Tuple{Int,Int}[]
    for i in 1:nx, j in (1, ny)
        (!dil[i, j] && !seen[i, j]) && (seen[i, j] = true; push!(stk, (i, j)))
    end
    for j in 1:ny, i in (1, nx)
        (!dil[i, j] && !seen[i, j]) && (seen[i, j] = true; push!(stk, (i, j)))
    end
    while !isempty(stk)
        (i, j) = pop!(stk)
        for (di, dj) in ((1, 0), (-1, 0), (0, 1), (0, -1))
            ii = i + di; jj = j + dj
            (1 <= ii <= nx && 1 <= jj <= ny && !dil[ii, jj] && !seen[ii, jj]) &&
                (seen[ii, jj] = true; push!(stk, (ii, jj)))
        end
    end
    fp = occ .| ((.!seen) .& (.!dil))
    return BitVector([fp[bx(xs[i]), by(ys[i])] for i in eachindex(xs)])
end

# ---- Enclosure reclass --------------------------------------------------------

# Exact interior: a background ("outside") point ENCLOSED by the cell (≥ `min_hits`
# of 8 rays hit cell tissue before the field edge) is inside the cell → reclassify
# `interior`. Removes residual interior background, including pockets reached
# through low-density channels. Mutates `class` in place; only `outside → interior`.
function _enclosure_fill!(class::Vector{String}, xs::Vector{Float64},
                          ys::Vector{Float64}; bin::Float64 = 0.2,
                          min_hits::Int = 6)
    N = length(class)
    x0, x1 = extrema(xs); y0, y1 = extrema(ys)
    nx = max(1, ceil(Int, (x1 - x0) / bin) + 1)
    ny = max(1, ceil(Int, (y1 - y0) / bin) + 1)
    bx(x) = clamp(floor(Int, (x - x0) / bin) + 1, 1, nx)
    by(y) = clamp(floor(Int, (y - y0) / bin) + 1, 1, ny)
    cell = falses(nx, ny)
    for i in 1:N
        (class[i] == "interior" || class[i] == "membrane") &&
            (cell[bx(xs[i]), by(ys[i])] = true)
    end
    dirs = ((1, 0), (-1, 0), (0, 1), (0, -1), (1, 1), (1, -1), (-1, 1), (-1, -1))
    for i in 1:N
        class[i] == "outside" || continue
        bi = bx(xs[i]); bj = by(ys[i]); hits = 0
        for (dx, dy) in dirs
            ci = bi; cj = bj
            while true
                ci += dx; cj += dy
                (1 <= ci <= nx && 1 <= cj <= ny) || break
                cell[ci, cj] && (hits += 1; break)
            end
        end
        hits >= min_hits && (class[i] = "interior")
    end
    return class
end

# ---- Gate + full path ---------------------------------------------------------

# Stages 1–2: KDE-valley density gate + footprint fill on the ORIGINAL cloud.
# Returns the per-point in-footprint mask (the tissue handed to the v1 geometry).
function _kde_valley_footprint(x::Vector{Float64}, y::Vector{Float64},
                               params::EdgeClassifyConfig)
    n = length(x)
    X = Matrix{Float64}(undef, 2, n)
    @inbounds for i in 1:n
        X[1, i] = x[i]; X[2, i] = y[i]
    end
    tree = NearestNeighbors.KDTree(X)
    sigma = params.KDE_SIGMA_NM / 1000          # nm → µm
    rho = _kde_density(X, tree, sigma; rmax_sigma = params.KDE_RMAX_SIGMA)
    rho_thr = 10^_kde_leftbase(log10.(rho .+ 1.0);
                               nbins = params.KDE_VALLEY_NBINS,
                               floorfrac = params.KDE_VALLEY_FLOORFRAC,
                               smooth = params.KDE_VALLEY_SMOOTH) - 1.0
    tissue = rho .>= rho_thr
    return _footprint_fill(x, y, tissue;
                           bin = params.FOOTPRINT_BIN_UM,
                           closing = params.FOOTPRINT_CLOSING_PX)
end

# Full kde_valley classification path. Called from `classify_emitters` when
# `params.METHOD == "kde_valley"`. Reproduces the validated genmab order
# gate(originals) → footprint → reflect → alpha → classify → enclosure.
function _classify_kde_valley(x::Vector{Float64}, y::Vector{Float64},
                              fov_um::NTuple{4,Float64},
                              params::EdgeClassifyConfig, t0::Float64;
                              out_dir, condition, cell,
                              write_artifacts::Bool, write_renders::Bool,
                              smld_input_meta)
    n = length(x)

    # Stages 1–2: KDE-valley gate + footprint on the original (pre-reflection) cloud.
    tfill = _kde_valley_footprint(x, y, params)
    any(tfill) || throw(ErrorException(
        "METHOD=\"kde_valley\": KDE-valley + footprint gate produced empty tissue " *
        "(KDE_SIGMA_NM=$(params.KDE_SIGMA_NM)); cloud too sparse or sigma too small"))

    # v1 outer-polygon geometry on the footprint SUBSET. Internal k-NN gate OFF
    # (RHO_K_THRESH = 0): the KDE gate already selected tissue. Reflection +
    # alpha-shape + point-in-polygon + membrane band happen inside on the subset.
    fidx = findall(tfill)
    subx = x[fidx]; suby = y[fidx]
    sub_params = EdgeClassifyConfig(
        METHOD            = _METHOD_OUTER_POLYGON,
        RHO_K_THRESH      = 0.0,
        K_LIST            = params.K_LIST,
        ALPHA_NM          = params.ALPHA_NM,
        REFLECT_RADIUS_NM = params.REFLECT_RADIUS_NM,
        MEMBRANE_NM       = params.MEMBRANE_NM,
        FOV_TRUNC_TOL_NM  = params.FOV_TRUNC_TOL_NM,
    )
    sub = classify_emitters(subx, suby; fov_um = fov_um, params = sub_params)

    # Scatter subset classification back to the full original array; off-footprint
    # points are background ("outside").
    class = fill("outside", n)
    inside_outer = falses(n)
    dist_to_outer = fill(NaN, n)
    @inbounds for (k, i) in enumerate(fidx)
        class[i] = sub.class[k]
        inside_outer[i] = sub.inside_outer[k]
        dist_to_outer[i] = sub.dist_to_outer_um[k]
    end

    # Stage 3: enclosure reclass on the full original cloud. Per topology contract
    # (i): folds enclosed background into `class == "interior"`; `inside_outer`/
    # `dist_to_outer_um` stay geometric (dist NaN for enclosure-promoted points).
    _enclosure_fill!(class, x, y;
                     bin = params.ENCLOSURE_BIN_UM,
                     min_hits = params.ENCLOSURE_MIN_HITS)

    # in_cell = topological membership: geometric-inside ∪ enclosure-recovered.
    in_cell = .!(class .== "outside")

    runtime_s = time() - t0
    result = EdgeClassificationResult(
        n, class, inside_outer, in_cell, dist_to_outer,
        sub.outer_polygon, sub.loops, sub.loop_diagnostics,
        params, fov_um, sub.truncated_sides, sub.n_reflected, runtime_s,
        nothing,
    )

    if write_artifacts
        leaf = joinpath(out_dir, condition, cell)
        mkpath(leaf)
        _write_artifacts(leaf, result;
                         condition = condition, cell = cell,
                         smld_input_meta = smld_input_meta,
                         x_um = x, y_um = y,
                         write_renders = write_renders)
    end

    return result
end

# ---- Validated-preset factory -------------------------------------------------

"""
    kde_valley_params(; sigma_nm=150.0, alpha_nm=600.0, reflect_radius_nm=1500.0,
                      membrane_nm=100.0, kwargs...) -> EdgeClassifyConfig

Construct an [`EdgeClassifyConfig`](@ref) with `METHOD = "kde_valley"` and the
**validated A431 dSTORM defaults** (σ = 150 nm, α = 600 nm, reflect = 1500 nm,
membrane = 100 nm). Prefer this over a raw `EdgeClassifyConfig(METHOD="kde_valley")`:
the validated `ALPHA_NM = 600` differs from the struct default (300), so a raw
constructor would silently use the wrong alpha. Extra `kwargs` override any field.

```julia
res = classify_emitters(x_um, y_um; fov_um, params = kde_valley_params())
```
"""
function kde_valley_params(; sigma_nm::Real = 150.0, alpha_nm::Real = 600.0,
                           reflect_radius_nm::Real = 1500.0,
                           membrane_nm::Real = 100.0, kwargs...)
    return EdgeClassifyConfig(; METHOD = _METHOD_KDE_VALLEY,
                              KDE_SIGMA_NM = Float64(sigma_nm),
                              ALPHA_NM = Float64(alpha_nm),
                              REFLECT_RADIUS_NM = Float64(reflect_radius_nm),
                              MEMBRANE_NM = Float64(membrane_nm),
                              kwargs...)
end
