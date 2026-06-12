"""
Internal: the two dispatched refinement hooks for `classify_emitters` and their
implementations.

- `_effective_polygon(loops, x, y, fov, cfg)` — `MaskCarveConfig` substitutes a
  density-mask carve of the v1 outer polygon (carve ⊆ v1) and returns its
  `MaskCarveDiagnostic` as `aux`.
- `_refine!(class, inside_outer, dist, x, y, fov, cfg)` — `GridHybridConfig`
  promotes near-boundary interior emitters to `:membrane`.
- `_enclosure_fill!(class, x, y, cfg)` — `KdeValleyConfig` enclosure pass (called
  directly from the kde_valley core, not a hook).

Validated math ported verbatim from `grid_hybrid.jl` / `mask_carve.jl` /
`edge_mask.jl`; only the parameter source + class representation change.
"""

# ===== mask_carve: effective-polygon hook =====================================

function _effective_polygon(loops, x, y, fov, cfg::MaskCarveConfig)
    return _build_mask_carve(loops[1], x, y, fov, cfg)
end

# ===== grid_hybrid: refinement hook ===========================================

function _refine!(class::Vector{Symbol}, inside_outer, dist::Vector{Float64},
                  x, y, fov, cfg::GridHybridConfig)
    grid_membrane = _grid_boundary_membrane_mask(x, y, fov, cfg)
    max_outer_dist_um = cfg.grid_outer_buffer_nm / 1000
    @inbounds for i in eachindex(class)
        if class[i] == :interior &&
           grid_membrane[i] &&
           isfinite(dist[i]) &&
           dist[i] <= max_outer_dist_um
            class[i] = :membrane
        end
    end
    return nothing
end

# ===== kde_valley: enclosure reclass ==========================================

# A background (:outside) point ENCLOSED by the cell (≥ `enclosure_min_hits` of 8
# rays hit cell tissue before the field edge) is folded into `:interior`. Mutates
# `class` in place; only `:outside → :interior`.
function _enclosure_fill!(class::Vector{Symbol}, xs::Vector{Float64},
                          ys::Vector{Float64}, cfg::KdeValleyConfig)
    bin = cfg.enclosure_bin_um
    min_hits = cfg.enclosure_min_hits
    N = length(class)
    x0, x1 = extrema(xs); y0, y1 = extrema(ys)
    nx = max(1, ceil(Int, (x1 - x0) / bin) + 1)
    ny = max(1, ceil(Int, (y1 - y0) / bin) + 1)
    bx(x) = clamp(floor(Int, (x - x0) / bin) + 1, 1, nx)
    by(y) = clamp(floor(Int, (y - y0) / bin) + 1, 1, ny)
    cell = falses(nx, ny)
    for i in 1:N
        (class[i] == :interior || class[i] == :membrane) &&
            (cell[bx(xs[i]), by(ys[i])] = true)
    end
    dirs = ((1, 0), (-1, 0), (0, 1), (0, -1), (1, 1), (1, -1), (-1, 1), (-1, -1))
    for i in 1:N
        class[i] == :outside || continue
        bi = bx(xs[i]); bj = by(ys[i]); hits = 0
        for (dx, dy) in dirs
            ci = bi; cj = bj
            while true
                ci += dx; cj += dy
                (1 <= ci <= nx && 1 <= cj <= ny) || break
                cell[ci, cj] && (hits += 1; break)
            end
        end
        hits >= min_hits && (class[i] = :interior)
    end
    return class
end

# ===== grid_hybrid helpers ====================================================

const _GRID_NBR4 = ((1, 0), (-1, 0), (0, 1), (0, -1))

function _gaussian_kernel(sigma_px::Float64)
    r = max(1, ceil(Int, 3 * sigma_px))
    offs = collect(-r:r)
    vals = exp.(-(offs .^ 2) ./ (2 * sigma_px^2))
    vals ./= sum(vals)
    return offs, vals
end

function _smooth_grid(A::Matrix{Float64}, sigma_px::Float64)
    offs, vals = _gaussian_kernel(sigma_px)
    ny, nx = size(A)
    B = zeros(Float64, ny, nx)
    C = zeros(Float64, ny, nx)
    @inbounds for y in 1:ny, x in 1:nx
        s = 0.0
        for (k, dx) in pairs(offs)
            xx = x + dx
            1 <= xx <= nx || continue
            s += vals[k] * A[y, xx]
        end
        B[y, x] = s
    end
    @inbounds for y in 1:ny, x in 1:nx
        s = 0.0
        for (k, dy) in pairs(offs)
            yy = y + dy
            1 <= yy <= ny || continue
            s += vals[k] * B[yy, x]
        end
        C[y, x] = s
    end
    return C
end

function _rasterize_points(x, y, fov::NTuple{4,Float64}, px_um::Float64)
    fxmin, fxmax, fymin, fymax = fov
    nx = max(1, ceil(Int, (fxmax - fxmin) / px_um))
    ny = max(1, ceil(Int, (fymax - fymin) / px_um))
    grid = zeros(Float64, ny, nx)
    @inbounds for i in eachindex(x)
        fxmin <= x[i] <= fxmax || continue
        fymin <= y[i] <= fymax || continue
        ix = clamp(floor(Int, (x[i] - fxmin) / px_um) + 1, 1, nx)
        iy = clamp(floor(Int, (y[i] - fymin) / px_um) + 1, 1, ny)
        grid[iy, ix] += 1.0
    end
    return grid, (xmin = fxmin, ymin = fymin, px_um = px_um, nx = nx, ny = ny)
end

function _keep_largest_component(mask::BitMatrix)
    ny, nx = size(mask)
    seen = falses(ny, nx)
    best = Tuple{Int,Int}[]
    @inbounds for y in 1:ny, x in 1:nx
        mask[y, x] || continue
        seen[y, x] && continue
        q = Tuple{Int,Int}[(y, x)]
        comp = Tuple{Int,Int}[]
        seen[y, x] = true
        head = 1
        while head <= length(q)
            cy, cx = q[head]
            head += 1
            push!(comp, (cy, cx))
            for (dy, dx) in _GRID_NBR4
                yy = cy + dy
                xx = cx + dx
                1 <= yy <= ny && 1 <= xx <= nx || continue
                mask[yy, xx] || continue
                seen[yy, xx] && continue
                seen[yy, xx] = true
                push!(q, (yy, xx))
            end
        end
        length(comp) > length(best) && (best = comp)
    end
    out = falses(ny, nx)
    for (y, x) in best
        out[y, x] = true
    end
    return out
end

function _fill_internal_holes(mask::BitMatrix)
    ny, nx = size(mask)
    exterior = falses(ny, nx)
    q = Tuple{Int,Int}[]
    for x in 1:nx
        if !mask[1, x]
            exterior[1, x] = true
            push!(q, (1, x))
        end
        if !mask[ny, x] && !exterior[ny, x]
            exterior[ny, x] = true
            push!(q, (ny, x))
        end
    end
    for y in 1:ny
        if !mask[y, 1] && !exterior[y, 1]
            exterior[y, 1] = true
            push!(q, (y, 1))
        end
        if !mask[y, nx] && !exterior[y, nx]
            exterior[y, nx] = true
            push!(q, (y, nx))
        end
    end
    head = 1
    while head <= length(q)
        cy, cx = q[head]
        head += 1
        for (dy, dx) in _GRID_NBR4
            yy = cy + dy
            xx = cx + dx
            1 <= yy <= ny && 1 <= xx <= nx || continue
            mask[yy, xx] && continue
            exterior[yy, xx] && continue
            exterior[yy, xx] = true
            push!(q, (yy, xx))
        end
    end
    return mask .| .!exterior
end

function _boundary_cells(mask::BitMatrix)
    ny, nx = size(mask)
    cells = Tuple{Int,Int}[]
    @inbounds for y in 1:ny, x in 1:nx
        mask[y, x] || continue
        is_boundary = y == 1 || y == ny || x == 1 || x == nx
        if !is_boundary
            for (dy, dx) in _GRID_NBR4
                if !mask[y + dy, x + dx]
                    is_boundary = true
                    break
                end
            end
        end
        is_boundary && push!(cells, (y, x))
    end
    return cells
end

function _grid_boundary_membrane_mask(x, y, fov::NTuple{4,Float64}, cfg::GridHybridConfig)
    px_um = cfg.grid_px_nm / 1000
    smooth_um = cfg.grid_smooth_nm / 1000
    membrane_um = cfg.membrane_nm / 1000
    counts, geom = _rasterize_points(x, y, fov, px_um)
    smooth = _smooth_grid(counts, smooth_um / px_um)
    nz = filter(>(0), vec(smooth))
    isempty(nz) && return falses(length(x))
    thr = max(quantile(nz, cfg.grid_mask_q),
              cfg.grid_mask_peak_frac * maximum(smooth))
    mask = BitMatrix(smooth .> thr)
    mask = _keep_largest_component(mask)
    mask = _fill_internal_holes(mask)
    boundary = _boundary_cells(mask)
    isempty(boundary) && return falses(length(x))

    B = Matrix{Float64}(undef, 2, length(boundary))
    @inbounds for (i, (yy, xx)) in enumerate(boundary)
        B[1, i] = geom.xmin + (xx - 0.5) * geom.px_um
        B[2, i] = geom.ymin + (yy - 0.5) * geom.px_um
    end
    tree = KDTree(B)
    out = falses(length(x))
    @inbounds for i in eachindex(x)
        ix = clamp(floor(Int, (x[i] - geom.xmin) / geom.px_um) + 1, 1, geom.nx)
        iy = clamp(floor(Int, (y[i] - geom.ymin) / geom.px_um) + 1, 1, geom.ny)
        mask[iy, ix] || continue
        _, ds = knn(tree, [x[i], y[i]], 1, true)
        ds[1] <= membrane_um && (out[i] = true)
    end
    return out
end

# ===== mask_carve helpers (self-contained, [nx, ny] orientation) ==============

function _mc_kde_grid(x_um::Vector{Float64}, y_um::Vector{Float64},
                     fov::NTuple{4,Float64}, pixel_um::Float64, sigma_um::Float64)
    fxmin, fxmax, fymin, fymax = fov
    nx = max(1, round(Int, (fxmax - fxmin) / pixel_um))
    ny = max(1, round(Int, (fymax - fymin) / pixel_um))
    g = zeros(Float64, nx, ny)
    @inbounds for i in eachindex(x_um)
        xi = x_um[i]; yi = y_um[i]
        (fxmin <= xi <= fxmax && fymin <= yi <= fymax) || continue
        cx = clamp(floor(Int, (xi - fxmin) / pixel_um) + 1, 1, nx)
        cy = clamp(floor(Int, (yi - fymin) / pixel_um) + 1, 1, ny)
        g[cx, cy] += 1.0
    end
    _mc_gaussian_blur_separable!(g, sigma_um / pixel_um)
    return g, nx, ny
end

function _mc_gaussian_blur_separable!(g::Matrix{Float64}, sigma_px::Float64)
    nx, ny = size(g)
    half = max(1, ceil(Int, 3 * sigma_px))
    kern = [exp(-0.5 * (k / sigma_px)^2) for k in -half:half]
    kern ./= sum(kern)
    tmp = similar(g)
    @inbounds for j in 1:ny, i in 1:nx
        s = 0.0
        for k in -half:half
            ii = clamp(i + k, 1, nx)
            s += kern[k + half + 1] * g[ii, j]
        end
        tmp[i, j] = s
    end
    @inbounds for j in 1:ny, i in 1:nx
        s = 0.0
        for k in -half:half
            jj = clamp(j + k, 1, ny)
            s += kern[k + half + 1] * tmp[i, jj]
        end
        g[i, j] = s
    end
    return g
end

function _mc_otsu_threshold(values::AbstractVector{Float64}; nbins::Int = 256)
    isempty(values) && return 0.0
    vmin, vmax = extrema(values)
    vmax <= vmin && return vmin
    counts = zeros(Int, nbins)
    @inbounds for v in values
        b = clamp(floor(Int, (v - vmin) / (vmax - vmin) * nbins) + 1, 1, nbins)
        counts[b] += 1
    end
    bin_centers = [vmin + (k - 0.5) * (vmax - vmin) / nbins for k in 1:nbins]
    total = sum(counts)
    sumT = sum(counts .* bin_centers)
    sumB = 0.0; wB = 0; max_var = -Inf; threshold = vmin
    @inbounds for k in 1:nbins
        wB += counts[k]
        wB == 0 && continue
        wF = total - wB
        wF == 0 && break
        sumB += counts[k] * bin_centers[k]
        mB = sumB / wB
        mF = (sumT - sumB) / wF
        between = wB * wF * (mB - mF)^2
        if between > max_var
            max_var = between
            threshold = vmin + k * (vmax - vmin) / nbins
        end
    end
    return threshold
end

function _mc_estimate_noise_floor_otsu(g::Matrix{Float64})
    nz = filter(x -> x > 0.0, vec(g))
    isempty(nz) && return 0.0
    log_nz = log10.(1 .+ 1e6 .* nz)
    τ_log = _mc_otsu_threshold(log_nz)
    τ_g = (10.0^τ_log - 1.0) / 1e6
    return max(τ_g, 1e-12)
end

function _mc_erode_3x3(mask::BitMatrix)
    nx, ny = size(mask)
    out = falses(nx, ny)
    @inbounds for j in 1:ny, i in 1:nx
        mask[i, j] || continue
        ok = true
        for dj in -1:1, di in -1:1
            ii = i + di; jj = j + dj
            if 1 <= ii <= nx && 1 <= jj <= ny
                mask[ii, jj] || (ok = false; break)
            else
                ok = false; break
            end
        end
        ok && (out[i, j] = true)
    end
    return out
end

function _mc_dilate_3x3(mask::BitMatrix)
    nx, ny = size(mask)
    out = falses(nx, ny)
    @inbounds for j in 1:ny, i in 1:nx
        mask[i, j] || continue
        for dj in -1:1, di in -1:1
            ii = i + di; jj = j + dj
            (1 <= ii <= nx && 1 <= jj <= ny) && (out[ii, jj] = true)
        end
    end
    return out
end

_mc_morph_open(mask::BitMatrix) = _mc_dilate_3x3(_mc_erode_3x3(mask))

function _mc_connected_components(mask::BitMatrix)
    nx, ny = size(mask)
    labels = zeros(Int, nx, ny)
    sizes = Int[]
    next_label = 0
    queue = Tuple{Int,Int}[]
    @inbounds for j in 1:ny, i in 1:nx
        if mask[i, j] && labels[i, j] == 0
            next_label += 1
            push!(sizes, 0)
            empty!(queue); push!(queue, (i, j))
            labels[i, j] = next_label
            while !isempty(queue)
                ci, cj = pop!(queue)
                sizes[next_label] += 1
                for (di, dj) in ((-1, 0), (1, 0), (0, -1), (0, 1))
                    ni = ci + di; nj = cj + dj
                    if 1 <= ni <= nx && 1 <= nj <= ny &&
                       mask[ni, nj] && labels[ni, nj] == 0
                        labels[ni, nj] = next_label
                        push!(queue, (ni, nj))
                    end
                end
            end
        end
    end
    return labels, sizes
end

function _mc_largest_cc(mask::BitMatrix)
    labels, sizes = _mc_connected_components(mask)
    isempty(sizes) && return falses(size(mask)...), 0
    best = argmax(sizes)
    out = falses(size(mask)...)
    nx, ny = size(mask)
    @inbounds for j in 1:ny, i in 1:nx
        labels[i, j] == best && (out[i, j] = true)
    end
    return out, sizes[best]
end

function _mc_drop_small_components(mask::BitMatrix, min_pix::Int)
    labels, sizes = _mc_connected_components(mask)
    isempty(sizes) && return falses(size(mask)...)
    nx, ny = size(mask)
    out = falses(nx, ny)
    @inbounds for j in 1:ny, i in 1:nx
        L = labels[i, j]
        L > 0 && sizes[L] >= min_pix && (out[i, j] = true)
    end
    return out
end

function _mc_fill_internal_holes_size(mask::BitMatrix, max_pix::Int)
    nx, ny = size(mask)
    inv = .!mask
    labels, sizes = _mc_connected_components(inv)
    touches_boundary = falses(length(sizes))
    @inbounds for j in 1:ny, i in 1:nx
        if (i == 1 || i == nx || j == 1 || j == ny) && labels[i, j] > 0
            touches_boundary[labels[i, j]] = true
        end
    end
    out = copy(mask)
    n_filled = 0; n_preserved = 0
    @inbounds for L in 1:length(sizes)
        touches_boundary[L] && continue
        if sizes[L] <= max_pix
            n_filled += 1
        else
            n_preserved += 1
        end
    end
    @inbounds for j in 1:ny, i in 1:nx
        L = labels[i, j]
        if L > 0 && !touches_boundary[L] && sizes[L] <= max_pix
            out[i, j] = true
        end
    end
    return out, n_filled, n_preserved
end

function _mc_find_top_left_cell(mask::BitMatrix)
    nx, ny = size(mask)
    @inbounds for j in 1:ny, i in 1:nx
        mask[i, j] && return (i, j)
    end
    return (0, 0)
end

function _mc_trace_outer_contour(mask::BitMatrix; max_iters_factor::Int = 10)
    nx, ny = size(mask)
    start = _mc_find_top_left_cell(mask)
    start == (0, 0) && return Tuple{Int,Int}[]
    moore = [(0, -1), (1, -1), (1, 0), (1, 1), (0, 1), (-1, 1), (-1, 0), (-1, -1)]
    contour = Tuple{Int,Int}[start]
    current = start
    backtrack = 7
    max_iters = max_iters_factor * (nx + ny)
    for _ in 1:max_iters
        found = false
        for k in 1:8
            d_idx = mod1(backtrack + k, 8)
            di, dj = moore[d_idx]
            ni, nj = current[1] + di, current[2] + dj
            if 1 <= ni <= nx && 1 <= nj <= ny && mask[ni, nj]
                current = (ni, nj)
                backtrack = mod1(d_idx + 4, 8)
                push!(contour, current)
                found = true
                break
            end
        end
        found || break
        if length(contour) >= 4 && current == start
            (contour[2] == contour[end-1] || length(contour) > 6) && break
        end
    end
    return contour
end

_mc_cell_to_xy(i, j, fov, pixel) = (fov[1] + (i - 0.5) * pixel,
                                    fov[3] + (j - 0.5) * pixel)

function _mc_polygon_to_mask(poly::AbstractVector{<:NTuple{2,Float64}},
                             fov::NTuple{4,Float64}, pixel_um::Float64,
                             nx::Int, ny::Int)
    fxmin, _, fymin, _ = fov
    mask = falses(nx, ny)
    @inbounds for j in 1:ny, i in 1:nx
        px = fxmin + (i - 0.5) * pixel_um
        py = fymin + (j - 0.5) * pixel_um
        mask[i, j] = _point_in_polygon(px, py, poly)
    end
    return mask
end

function _dist_to_polygon_polyline(qx::Float64, qy::Float64,
                                   poly::AbstractVector{<:NTuple{2,Float64}})
    n = length(poly); best = Inf
    @inbounds for i in 1:n
        a = poly[i]; b = poly[i == n ? 1 : i + 1]
        dx = b[1] - a[1]; dy = b[2] - a[2]
        L2 = dx*dx + dy*dy
        if L2 < 1e-18
            d = hypot(qx - a[1], qy - a[2])
            d < best && (best = d); continue
        end
        t = clamp(((qx - a[1]) * dx + (qy - a[2]) * dy) / L2, 0.0, 1.0)
        cx = a[1] + t*dx; cy = a[2] + t*dy
        d = hypot(qx - cx, qy - cy)
        d < best && (best = d)
    end
    return best
end

function _mc_symmetric_dist_summary(p1::AbstractVector{<:NTuple{2,Float64}},
                                    p2::AbstractVector{<:NTuple{2,Float64}})
    if isempty(p1) || isempty(p2)
        return (NaN, NaN)
    end
    a = [_dist_to_polygon_polyline(v[1], v[2], p2) for v in p1]
    b = [_dist_to_polygon_polyline(v[1], v[2], p1) for v in p2]
    all_d = vcat(a, b)
    return (median(all_d), quantile(all_d, 0.95))
end

function _polygon_area_abs(poly::AbstractVector{<:NTuple{2,Float64}})
    n = length(poly)
    n < 3 && return 0.0
    s = 0.0
    @inbounds for i in 1:n
        x1, y1 = poly[i]
        x2, y2 = poly[mod1(i + 1, n)]
        s += x1*y2 - x2*y1
    end
    return abs(s) / 2
end

# Returns (effective_polygon, MaskCarveDiagnostic). On any degeneracy returns
# (v1_polygon, diagnostic with applied=false).
function _build_mask_carve(v1_polygon::AbstractVector{<:NTuple{2,Float64}},
                           x_um::Vector{Float64}, y_um::Vector{Float64},
                           fov_um::NTuple{4,Float64}, cfg::MaskCarveConfig)
    σ       = cfg.sigma_um
    knoise  = cfg.k_noise
    pixel   = cfg.pixel_um
    minfrac = cfg.min_component_frac
    fill_um2 = cfg.fill_hole_max_um2

    px_area = pixel * pixel
    fill_max_pix = max(1, round(Int, fill_um2 / px_area))

    v1_area_um2 = _polygon_area_abs(v1_polygon)
    fallback_diag(reason::String) = MaskCarveDiagnostic(
        false, reason, σ, knoise, pixel, minfrac, fill_um2,
        v1_area_um2, v1_area_um2, 0.0, 0.0, 0.0, 0.0, 0.0,
        0, 0, length(v1_polygon))

    g, nx, ny = _mc_kde_grid(x_um, y_um, fov_um, pixel, σ)
    pmax = maximum(g)
    pmax <= 0 && return v1_polygon, fallback_diag("empty_density_grid")
    noise_floor = _mc_estimate_noise_floor_otsu(g)
    τ = max(knoise * noise_floor, 1e-6 * pmax)
    d_mask = falses(nx, ny)
    @inbounds for j in 1:ny, i in 1:nx
        d_mask[i, j] = g[i, j] >= τ
    end
    any(d_mask) || return v1_polygon, fallback_diag("empty_d_mask")

    d_mask = _mc_morph_open(d_mask)
    any(d_mask) || return v1_polygon, fallback_diag("empty_after_morph_open")

    _, sizes = _mc_connected_components(d_mask)
    isempty(sizes) && return v1_polygon, fallback_diag("no_components")
    largest = maximum(sizes)
    min_pix = max(8, round(Int, minfrac * largest))
    d_mask = _mc_drop_small_components(d_mask, min_pix)

    d_filled, n_filled, n_preserved = _mc_fill_internal_holes_size(d_mask, fill_max_pix)

    v1_mask = _mc_polygon_to_mask(v1_polygon, fov_um, pixel, nx, ny)
    carve_mask = falses(nx, ny)
    n_carve = 0
    @inbounds for j in 1:ny, i in 1:nx
        if v1_mask[i, j] && d_filled[i, j]
            carve_mask[i, j] = true
            n_carve += 1
        end
    end
    n_carve == 0 && return v1_polygon, fallback_diag("empty_intersection")

    cc, _ = _mc_largest_cc(carve_mask)
    contour = _mc_trace_outer_contour(cc)
    length(contour) < 4 && return v1_polygon, fallback_diag("trace_too_short")
    contour[end] == contour[1] || return v1_polygon, fallback_diag("trace_not_closed")
    poly = NTuple{2,Float64}[_mc_cell_to_xy(i, j, fov_um, pixel) for (i, j) in contour]

    v1_area = _polygon_area_abs(v1_polygon)
    carve_area = _polygon_area_abs(poly)
    v_only_pix = 0; c_only_pix = 0
    carve_rast = _mc_polygon_to_mask(poly, fov_um, pixel, nx, ny)
    @inbounds for j in 1:ny, i in 1:nx
        iv = v1_mask[i, j]; ic = carve_rast[i, j]
        if iv && !ic; v_only_pix += 1
        elseif ic && !iv; c_only_pix += 1; end
    end
    v_only_a = v_only_pix * px_area
    c_only_a = c_only_pix * px_area
    med_d, p95_d = _mc_symmetric_dist_summary(v1_polygon, poly)

    diag = MaskCarveDiagnostic(
        true, "", σ, knoise, pixel, minfrac, fill_um2,
        v1_area, carve_area, carve_area - v1_area,
        v_only_a, c_only_a, med_d, p95_d,
        n_filled, n_preserved, length(poly))
    return poly, diag
end
