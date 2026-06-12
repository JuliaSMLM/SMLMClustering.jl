"""
Internal: density gates. The multi-K k-NN tissue mask lives in `geometry.jl`
(`_tissue_mask`); this file holds the KDE-valley gate used by `KdeValleyConfig`.

Validated functions ported verbatim from the genmab production wrapper
(`paper-genmab-hexabody/src/edge_mask.jl`); only the parameter source changes
(config fields instead of a flat dict).
"""

# Continuous Gaussian-KDE density at each point (fixed bandwidth `sigma`, µm).
# Per-point range query: memory-safe on dense clouds, unbiased at edges/clusters
# (unlike k-NN ρ_K). Self-contribution subtracted.
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

# Left base of the dominant (cell) mode: scan left from the global peak of the
# smoothed log-density histogram to the first bin below `floorfrac` of the peak.
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

# Binary dilate (`dilate=true`) / erode (`dilate=false`) by a square element of radius `r`.
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

# Footprint of kept tissue: rasterize → dilate-seal thin necks → flood-fill
# enclosed holes from the grid border. Returns a per-point in-footprint mask.
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

# Stages 1–2 of kde_valley: KDE-valley density gate + footprint fill on the
# ORIGINAL cloud. Returns the per-point in-footprint mask handed to the polygon core.
function _kde_valley_footprint(x::Vector{Float64}, y::Vector{Float64},
                               cfg::KdeValleyConfig)
    n = length(x)
    X = Matrix{Float64}(undef, 2, n)
    @inbounds for i in 1:n
        X[1, i] = x[i]; X[2, i] = y[i]
    end
    tree = NearestNeighbors.KDTree(X)
    sigma = cfg.sigma_nm / 1000          # nm → µm
    rho = _kde_density(X, tree, sigma; rmax_sigma = cfg.rmax_sigma)
    rho_thr = 10^_kde_leftbase(log10.(rho .+ 1.0);
                               nbins = cfg.valley_nbins,
                               floorfrac = cfg.valley_floorfrac,
                               smooth = cfg.valley_smooth) - 1.0
    tissue = rho .>= rho_thr
    return _footprint_fill(x, y, tissue;
                           bin = cfg.footprint_bin_um,
                           closing = cfg.footprint_closing_px)
end
