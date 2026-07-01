"""
Internal: compute per-loop diagnostics for `loop_diagnostics.csv`.

Density at vertices is evaluated against the **originals-only** KDTree —
reflected/augmented points are used only for polygon construction, never
for diagnostics.
"""

const _DIAG_DENSITY_K = 128  # K used for vertex-density diagnostic
const _DIAG_DENSE_THRESH_FRAC = 0.8  # provisional heuristic_type threshold

function _heuristic_type(used_in_outer::Bool, frac_in_fov::Float64,
                         frac_dense::Float64)::String
    used_in_outer && return "outer"                 # a ring of the published mask (any cell)
    frac_in_fov < 1.0 && return "fov_crossing"
    frac_dense >= _DIAG_DENSE_THRESH_FRAC ? "interior_dense" : "interior_sparse"
end

function _compute_loop_diagnostics(
    loops::Vector{Vector{NTuple{2,Float64}}},
    x_um::Vector{Float64}, y_um::Vector{Float64},
    Xorig::Matrix{Float64},
    fov_um::NTuple{4,Float64},
    rho_thresh::Float64,
    cells,
)
    fxmin, fxmax, fymin, fymax = fov_um
    n_orig = length(x_um)
    # Rings the published multi-cell mask actually uses (each cell's outer + any holes),
    # as vertex sets — used below to mark `used_in_outer` independent of start vertex /
    # orientation. A loop that build_mask split (or dropped as sub-cutoff debris) won't
    # match, by design.
    mask_rings = Set{Set{NTuple{2,Float64}}}()
    for c in cells
        push!(mask_rings, Set(c.outer))
        for h in c.holes
            push!(mask_rings, Set(h))
        end
    end
    tree = NearestNeighbors.KDTree(Xorig)
    K = min(_DIAG_DENSITY_K, max(n_orig, 1))   # clamp so small clouds don't throw on knn
    inv_pi = 1 / π

    diags = Vector{LoopDiagnostic}(undef, length(loops))
    for (lid, verts) in enumerate(loops)
        nv = length(verts)
        area = abs(_polygon_area(verts))

        # n_emitters_inside (originals strictly inside)
        n_inside = 0
        @inbounds for i in 1:n_orig
            _point_in_polygon(x_um[i], y_um[i], verts) && (n_inside += 1)
        end

        # vertex density vs originals KDTree
        rhos = Vector{Float64}(undef, nv)
        in_fov_count = 0
        @inbounds for vi in 1:nv
            vx, vy = verts[vi]
            (fxmin <= vx <= fxmax && fymin <= vy <= fymax) && (in_fov_count += 1)
            _, dists = NearestNeighbors.knn(tree, [vx, vy], K, true)
            d = dists[end]
            rhos[vi] = (K - 1) * inv_pi / (d * d)
        end
        frac_dense = count(>=(rho_thresh), rhos) / nv
        med_rho = Statistics.median(rhos)
        frac_in_fov = in_fov_count / nv

        # used_in_outer: this loop became a ring of the published multi-cell mask (a
        # cell outer or hole). Loops dropped as sub-cutoff debris (or split by
        # build_mask) are diagnostic-only. Relies on float-EXACT vertex equality —
        # holds because build_mask reuses the loop's untransformed vertices; revisit
        # if build_mask ever resamples/simplifies rings.
        used_in_outer = Set(verts) in mask_rings
        diags[lid] = LoopDiagnostic(
            lid, nv, area, n_inside, frac_in_fov, frac_dense, med_rho,
            used_in_outer,
            _heuristic_type(used_in_outer, frac_in_fov, frac_dense),
        )
    end
    return diags
end
