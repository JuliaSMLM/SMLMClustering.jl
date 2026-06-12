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

function _grid_boundary_membrane_mask(x, y, fov::NTuple{4,Float64},
                                      params::EdgeClassifyConfig)
    px_um = params.GRID_PX_NM / 1000
    smooth_um = params.GRID_SMOOTH_NM / 1000
    membrane_um = params.MEMBRANE_NM / 1000
    counts, geom = _rasterize_points(x, y, fov, px_um)
    smooth = _smooth_grid(counts, smooth_um / px_um)
    nz = filter(>(0), vec(smooth))
    isempty(nz) && return falses(length(x))
    thr = max(quantile(nz, params.GRID_MASK_Q),
              params.GRID_MASK_PEAK_FRAC * maximum(smooth))
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

function _apply_grid_hybrid!(class::Vector{String}, x, y,
                             dist_to_outer::Vector{Float64},
                             fov::NTuple{4,Float64},
                             params::EdgeClassifyConfig)
    grid_membrane = _grid_boundary_membrane_mask(x, y, fov, params)
    max_outer_dist_um = params.GRID_OUTER_BUFFER_NM / 1000
    @inbounds for i in eachindex(class)
        if class[i] == "interior" &&
           grid_membrane[i] &&
           isfinite(dist_to_outer[i]) &&
           dist_to_outer[i] <= max_outer_dist_um
            class[i] = "membrane"
        end
    end
    return class
end
