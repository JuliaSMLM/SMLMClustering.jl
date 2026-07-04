"""
Reusable polygon-region types for SMLM masks: a cell footprint as an outer ring
plus optional internal holes, and a multi-cell mask as a vector of those. Shared
by `EdgeClassify` (the published mask) and `Hopkins` (the observation window).
Pure Julia — point-in-region reuses `_point_in_polygon` (utils.jl); no geometry
dependencies. A `GeometryBasics` bridge is provided as an optional package
extension (loaded only if the user has GeometryBasics).
"""

"""
    CellPolygon(outer, holes = [])

One cell's footprint: an `outer` boundary ring (a `Vector` of `(x, y)` µm tuples,
stored open / not repeating the first vertex) plus optional internal `holes`
(each a ring). `holes` is empty unless the mask was built with
`keep_internal = true`.
"""
struct CellPolygon
    outer::Vector{NTuple{2,Float64}}
    holes::Vector{Vector{NTuple{2,Float64}}}
end
CellPolygon(outer::Vector{NTuple{2,Float64}}) =
    CellPolygon(outer, Vector{NTuple{2,Float64}}[])

"""
    MultiCellMask

`Vector{CellPolygon}` — the edge mask of a field of view: one `CellPolygon` per
distinct cell. Returned by `classify_emitters` (`info.cells` /
`metadata["edge_cells"]`) and accepted as a Hopkins observation window. Cells are
ordered largest-first, so `mask[1]` is the dominant cell.
"""
const MultiCellMask = Vector{CellPolygon}

# Shoelace area magnitude of a ring.
function _ring_area(ring::AbstractVector{NTuple{2,Float64}})
    n = length(ring)
    n < 3 && return 0.0
    s = 0.0
    @inbounds for i in 1:n
        j = i == n ? 1 : i + 1
        s += ring[i][1] * ring[j][2] - ring[j][1] * ring[i][2]
    end
    return abs(s) / 2
end

# Split a closed vertex ring into SIMPLE sub-rings at any repeated vertex (a
# self-touching figure-8 boundary → its separate lobes). A repeated coordinate is
# a pinch point; each maximal repeat-free span of ≥ 3 vertices becomes a sub-ring.
function _split_simple(loop::AbstractVector{NTuple{2,Float64}})
    subs = Vector{Vector{NTuple{2,Float64}}}()
    path = NTuple{2,Float64}[]
    pos = Dict{NTuple{2,Float64},Int}()
    @inbounds for v in loop
        idx = get(pos, v, 0)
        if idx != 0
            sub = path[idx:end]
            length(sub) >= 3 && push!(subs, copy(sub))
            for k in (idx + 1):length(path)
                delete!(pos, path[k])
            end
            resize!(path, idx)
        else
            push!(path, v)
            pos[v] = length(path)
        end
    end
    length(path) >= 3 && push!(subs, copy(path))
    return subs
end

# Robust nesting test: ring `a` is inside ring `b` when the MAJORITY of a's
# vertices fall inside b. A majority vote (not a single vertex) tolerates the
# shared pinch vertices that sit on b's boundary after a figure-8 split.
function _ring_in_ring(a::AbstractVector{NTuple{2,Float64}},
                       b::AbstractVector{NTuple{2,Float64}})
    inside = 0
    @inbounds for v in a
        _point_in_polygon(v[1], v[2], b) && (inside += 1)
    end
    return inside > length(a) ÷ 2
end

"""
    build_mask(loops; keep_internal = false, min_cell_frac = 1/3, min_hole_frac = 0) -> MultiCellMask

Assemble a [`MultiCellMask`](@ref) from raw alpha-shape boundary `loops` (closed
vertex rings):

1. **Split** each loop into simple sub-rings at self-touch / shared vertices.
2. **Nest** by depth (how many other rings contain it): even depth ⇒ a cell
   outer, odd depth ⇒ an internal region (a real void *or* a pinch-enclosed
   region).
3. **Group** each cell outer with the internal regions directly inside it — only
   when `keep_internal = true`; otherwise cells are solid (`holes` empty). With
   `min_hole_frac > 0`, holes smaller than `min_hole_frac × (that cell's outer area)`
   are dropped (filled back into the interior), so only real voids above that scale
   survive — sub-cell texture such as inter-cluster gaps no longer fragments the cell.
4. **Drop debris**: remove any cell whose outer area is below
   `min_cell_frac × (largest cell area)`; `min_cell_frac = 0` keeps everything.

Cells are returned largest-first.
"""
function build_mask(loops::AbstractVector; keep_internal::Bool = false,
                    min_cell_frac::Real = 1//3, min_hole_frac::Real = 0)
    (0 <= min_cell_frac < 1) ||
        throw(ArgumentError("min_cell_frac must be in [0,1); got $min_cell_frac"))
    (0 <= min_hole_frac < 1) ||
        throw(ArgumentError("min_hole_frac must be in [0,1); got $min_hole_frac"))
    simple = Vector{Vector{NTuple{2,Float64}}}()
    for L in loops
        for s in _split_simple(L)
            # Public-API guards (the internal pipeline never emits these): reject
            # non-finite coordinates, and drop degenerate rings (< 3 vertices or zero
            # signed area) that cannot bound a cell.
            length(s) >= 3 || continue
            all(p -> isfinite(p[1]) && isfinite(p[2]), s) ||
                throw(ArgumentError("build_mask: ring has non-finite coordinates"))
            abs(_ring_area(s)) > 0 || continue
            push!(simple, s)
        end
    end
    isempty(simple) && return CellPolygon[]
    n = length(simple)
    depth = zeros(Int, n)
    for i in 1:n
        @inbounds for j in 1:n
            i == j && continue
            _ring_in_ring(simple[i], simple[j]) && (depth[i] += 1)
        end
    end
    cells = CellPolygon[]
    for i in 1:n
        iseven(depth[i]) || continue
        hs = Vector{Vector{NTuple{2,Float64}}}()
        if keep_internal
            # Drop holes below min_hole_frac × this cell's outer area: separates
            # real internal voids (kept) from sub-cell texture like inter-cluster
            # gaps carved by the adaptive α-shape (filled → absorbed into interior).
            hole_thr = min_hole_frac > 0 ? min_hole_frac * _ring_area(simple[i]) : 0.0
            for j in 1:n
                depth[j] == depth[i] + 1 || continue
                _ring_in_ring(simple[j], simple[i]) || continue
                _ring_area(simple[j]) >= hole_thr && push!(hs, simple[j])
            end
        end
        push!(cells, CellPolygon(simple[i], hs))
    end
    isempty(cells) && return cells
    maxa = maximum(_ring_area(c.outer) for c in cells)
    if min_cell_frac > 0 && maxa > 0
        thr = min_cell_frac * maxa
        cells = [c for c in cells if _ring_area(c.outer) >= thr]
    end
    sort!(cells; by = c -> -_ring_area(c.outer))
    return cells
end

"""
    in_region(x, y, mask::MultiCellMask) -> Bool

True when `(x, y)` lies in any cell of `mask`: inside a cell's outer ring and not
inside one of that cell's holes.
"""
function in_region(x::Real, y::Real, mask::AbstractVector{CellPolygon})
    xf = Float64(x); yf = Float64(y)
    @inbounds for cell in mask
        _point_in_polygon(xf, yf, cell.outer) || continue
        inhole = false
        for h in cell.holes
            if _point_in_polygon(xf, yf, h)
                inhole = true
                break
            end
        end
        inhole || return true
    end
    return false
end

"""
    region_area(mask::MultiCellMask) -> Float64

Total occupied area in µm²: Σ cell outer areas − Σ hole areas.
"""
function region_area(mask::AbstractVector{CellPolygon})
    a = 0.0
    @inbounds for cell in mask
        a += _ring_area(cell.outer)
        for h in cell.holes
            a -= _ring_area(h)
        end
    end
    return a
end
