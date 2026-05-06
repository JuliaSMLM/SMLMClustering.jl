"""
Internal: mask_carve method (provisional/opt-in). Builds an effective
outer polygon by carving v1 inward using a heterogeneity-robust density
mask. Never expands v1 outward.

Pipeline (self-contained — does not share matrices with grid_hybrid; uses
[nx, ny] orientation throughout, where g[i, j] is x = fxmin + (i-0.5)*pixel,
y = fymin + (j-0.5)*pixel):

    1. KDE: rasterize emitters at MASK_CARVE_PIXEL_UM, separable Gaussian
       blur at MASK_CARVE_SIGMA_UM.
    2. Otsu noise floor: 2-class split on log10(1 + 1e6·g) over positive
       pixels only — robust regardless of cell-area fraction.
    3. Threshold: g ≥ MASK_CARVE_K_NOISE × noise_floor.
    4. Morphological open (3×3 erode then dilate).
    5. Drop CCs smaller than MASK_CARVE_MIN_COMPONENT_FRAC × largest.
    6. Fill internal holes whose area ≤ MASK_CARVE_FILL_HOLE_MAX_UM2
       (preserves legitimate large voids).
    7. Rasterize v1_polygon on same grid; intersect.
    8. Largest CC of intersection → Moore-neighbor 8-conn outer-boundary
       walk → polygon.
    9. Fallback: if any stage yields an empty/degenerate result, return
       v1_polygon and set diagnostic.applied = false with a reason.

Carve_polygon ⊆ v1_polygon by construction (any rounding-introduced
violations should be vanishingly small; reported via
`carve_only_area_um2`).
"""

# ---------- Grid / KDE helpers (self-contained, [nx, ny] orientation) ----

function _mc_kde_grid(x_um::Vector{Float64}, y_um::Vector{Float64},
                     fov::NTuple{4,Float64}, pixel_um::Float64, sigma_um::Float64)
    fxmin, fxmax, fymin, fymax = fov
    nx = max(1, round(Int, (fxmax-fxmin)/pixel_um))
    ny = max(1, round(Int, (fymax-fymin)/pixel_um))
    g = zeros(Float64, nx, ny)
    @inbounds for i in eachindex(x_um)
        xi = x_um[i]; yi = y_um[i]
        (fxmin <= xi <= fxmax && fymin <= yi <= fymax) || continue
        cx = clamp(floor(Int, (xi-fxmin)/pixel_um)+1, 1, nx)
        cy = clamp(floor(Int, (yi-fymin)/pixel_um)+1, 1, ny)
        g[cx, cy] += 1.0
    end
    _mc_gaussian_blur_separable!(g, sigma_um/pixel_um)
    return g, nx, ny
end

function _mc_gaussian_blur_separable!(g::Matrix{Float64}, sigma_px::Float64)
    nx, ny = size(g)
    half = max(1, ceil(Int, 3*sigma_px))
    kern = [exp(-0.5*(k/sigma_px)^2) for k in -half:half]
    kern ./= sum(kern)
    tmp = similar(g)
    @inbounds for j in 1:ny, i in 1:nx
        s = 0.0
        for k in -half:half
            ii = clamp(i+k, 1, nx)
            s += kern[k+half+1] * g[ii, j]
        end
        tmp[i, j] = s
    end
    @inbounds for j in 1:ny, i in 1:nx
        s = 0.0
        for k in -half:half
            jj = clamp(j+k, 1, ny)
            s += kern[k+half+1] * tmp[i, jj]
        end
        g[i, j] = s
    end
    return g
end

# ---------- Otsu --------------------------------------------------------

function _mc_otsu_threshold(values::AbstractVector{Float64}; nbins::Int=256)
    isempty(values) && return 0.0
    vmin, vmax = extrema(values)
    vmax <= vmin && return vmin
    counts = zeros(Int, nbins)
    @inbounds for v in values
        b = clamp(floor(Int, (v - vmin)/(vmax - vmin) * nbins) + 1, 1, nbins)
        counts[b] += 1
    end
    bin_centers = [vmin + (k - 0.5)*(vmax - vmin)/nbins for k in 1:nbins]
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
            # threshold = right edge of this bin
            threshold = vmin + k*(vmax - vmin)/nbins
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

# ---------- Morphology ([nx, ny] BitMatrix) -----------------------------

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

# ---------- Connected components ([nx, ny], 4-connectivity) -------------

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
                for (di, dj) in ((-1,0),(1,0),(0,-1),(0,1))
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

# Fill internal holes (non-FOV-boundary-touching) up to max_pix in size.
# Returns (filled_mask, n_holes_filled, n_holes_preserved).
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

# ---------- Polygonize: Moore-neighbor 8-conn outer-boundary walk -------

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
    moore = [(0,-1), (1,-1), (1,0), (1,1), (0,1), (-1,1), (-1,0), (-1,-1)]
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

_mc_cell_to_xy(i, j, fov, pixel) = (fov[1] + (i - 0.5)*pixel,
                                     fov[3] + (j - 0.5)*pixel)

# ---------- Polygon → mask rasterization on the same grid ----------------

function _mc_polygon_to_mask(poly::AbstractVector{<:NTuple{2,Float64}},
                              fov::NTuple{4,Float64}, pixel_um::Float64,
                              nx::Int, ny::Int)
    fxmin, _, fymin, _ = fov
    mask = falses(nx, ny)
    @inbounds for j in 1:ny, i in 1:nx
        px = fxmin + (i - 0.5)*pixel_um
        py = fymin + (j - 0.5)*pixel_um
        mask[i, j] = _point_in_polygon(px, py, poly)
    end
    return mask
end

# ---------- Per-vertex polygon-to-polygon distance (uses _dist_to_polygon) ----

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

# Distance from point to polygon polyline (closed). Independent of inside/outside.
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
        t = clamp(((qx - a[1])*dx + (qy - a[2])*dy) / L2, 0.0, 1.0)
        cx = a[1] + t*dx; cy = a[2] + t*dy
        d = hypot(qx - cx, qy - cy)
        d < best && (best = d)
    end
    return best
end

# ---------- Top-level builder -------------------------------------------

"""
    _build_mask_carve(v1_polygon, x_um, y_um, fov_um, params)
        -> (effective_polygon, MaskCarveDiagnostic)

Returns the carve polygon that should replace the v1 outer polygon for
classification, plus a diagnostic record. On any degeneracy the function
returns `v1_polygon` itself with `applied = false` and a reason — the
caller can use the result without re-checking.
"""
function _build_mask_carve(v1_polygon::AbstractVector{<:NTuple{2,Float64}},
                            x_um::Vector{Float64}, y_um::Vector{Float64},
                            fov_um::NTuple{4,Float64},
                            params::EdgeClassifyParams)
    σ      = params.MASK_CARVE_SIGMA_UM
    knoise = params.MASK_CARVE_K_NOISE
    pixel  = params.MASK_CARVE_PIXEL_UM
    minfrac = params.MASK_CARVE_MIN_COMPONENT_FRAC
    fill_um2 = params.MASK_CARVE_FILL_HOLE_MAX_UM2

    px_area = pixel * pixel
    fill_max_pix = max(1, round(Int, fill_um2 / px_area))

    # Fallback diagnostic: the effective polygon IS v1 (no carve happened),
    # so report carve_area = v1_area, all deltas/distances = 0, and
    # n_carve_polygon_pts = length(v1_polygon). Keeps the diagnostic
    # internally consistent with the result's effective polygon.
    v1_area_um2 = _polygon_area_abs(v1_polygon)
    fallback_diag(reason::String) = MaskCarveDiagnostic(
        false, reason, σ, knoise, pixel, minfrac, fill_um2,
        v1_area_um2, v1_area_um2, 0.0, 0.0, 0.0, 0.0, 0.0,
        0, 0, length(v1_polygon))

    # 1-3. KDE → noise floor → threshold.
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

    # 4. Morph open.
    d_mask = _mc_morph_open(d_mask)
    any(d_mask) || return v1_polygon, fallback_diag("empty_after_morph_open")

    # 5. Drop small CCs.
    _, sizes = _mc_connected_components(d_mask)
    isempty(sizes) && return v1_polygon, fallback_diag("no_components")
    largest = maximum(sizes)
    min_pix = max(8, round(Int, minfrac * largest))
    d_mask = _mc_drop_small_components(d_mask, min_pix)

    # 6. Fill internal holes (size-thresholded).
    d_filled, n_filled, n_preserved = _mc_fill_internal_holes_size(d_mask, fill_max_pix)

    # 7. v1 mask + intersection.
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

    # 8. Largest CC + polygonize.
    cc, _ = _mc_largest_cc(carve_mask)
    contour = _mc_trace_outer_contour(cc)
    length(contour) < 4 && return v1_polygon, fallback_diag("trace_too_short")
    # Closure invariant: a stable Moore-boundary walk ends with contour[end]
    # == contour[1] (the start cell is re-visited and pushed). If the walk
    # exited on max-iterations without closing, the trace is open — never
    # use it as an effective polygon.
    contour[end] == contour[1] || return v1_polygon, fallback_diag("trace_not_closed")
    poly = NTuple{2,Float64}[_mc_cell_to_xy(i, j, fov_um, pixel) for (i, j) in contour]

    # Areas + distances for diagnostic.
    v1_area = _polygon_area_abs(v1_polygon)
    carve_area = _polygon_area_abs(poly)
    # Direct rasterized v1\carve and carve\v1 (carve_only should ≈ 0)
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

# Helper: |signed polygon area| in µm². Co-located here to avoid touching
# geometry.jl's _polygon_area which returns signed.
function _polygon_area_abs(poly::AbstractVector{<:NTuple{2,Float64}})
    n = length(poly)
    n < 3 && return 0.0
    s = 0.0
    @inbounds for i in 1:n
        x1, y1 = poly[i]
        x2, y2 = poly[mod1(i+1, n)]
        s += x1*y2 - x2*y1
    end
    return abs(s) / 2
end
