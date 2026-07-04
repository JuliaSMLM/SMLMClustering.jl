"""
    classify_emitters(smld::BasicSMLD, cfg::AbstractEdgeClassifyConfig) -> (smld, info)
    classify_emitters(x_um, y_um, cfg::AbstractEdgeClassifyConfig; fov_um) -> info

Classify each emitter as `:outside`, `:membrane`, or `:interior`. The concrete
config type selects the **cell gate** by dispatch (`OuterPolygonConfig`,
`KdeValleyConfig`); the rest of the pipeline is shared.

The SMLD method follows the package `(out, Info)` convention: it returns the smld
with the published multi-cell mask threaded into `smld.metadata["edge_cells"]` (and
the dominant cell's outer ring into `smld.metadata["edge_outer_polygon"]` for
back-compat / Hopkins `region=:metadata`), plus an [`EdgeClassifyInfo`](@ref). The
per-emitter class lives **only** in `info.class` (with the `in_cell` / `interior_mask`
accessors) — it is deliberately *not* mirrored into the metadata, because a
per-emitter side-list desyncs the moment a downstream step subsets emitters. The
coordinate method is the computational core and returns the `info` directly;
`fov_um = (xmin, xmax, ymin, ymax)` in µm.
"""
function classify_emitters(smld::SMLMData.BasicSMLD, cfg::AbstractEdgeClassifyConfig)
    n = length(smld.emitters)
    x = Vector{Float64}(undef, n); y = Vector{Float64}(undef, n)
    @inbounds for i in 1:n
        x[i] = smld.emitters[i].x; y[i] = smld.emitters[i].y
    end
    fov = (Float64(smld.camera.pixel_edges_x[1]), Float64(smld.camera.pixel_edges_x[end]),
           Float64(smld.camera.pixel_edges_y[1]), Float64(smld.camera.pixel_edges_y[end]))
    info = classify_emitters(x, y, cfg; fov_um = fov)
    meta = copy(smld.metadata)
    # Only GEOMETRY is mirrored (safe under emitter-subsetting). The per-emitter class
    # is NOT — it lives in info.class (a per-emitter side-list in metadata would desync
    # the moment a downstream step filters/subsets emitters). Delete any stale class key
    # (from a prior classify pass on this SMLD) so it can't linger and silently desync.
    delete!(meta, "edge_classify_class")
    meta["edge_cells"] = info.cells                  # multi-cell mask (primary downstream output)
    meta["edge_outer_polygon"] = info.outer_polygon  # dominant cell outer (back-compat / Hopkins)
    smld_out = SMLMData.BasicSMLD(smld.emitters, smld.camera, smld.n_frames,
                                  smld.n_datasets, meta)
    return smld_out, info
end

function classify_emitters(x_um::AbstractVector{<:Real}, y_um::AbstractVector{<:Real},
                           cfg::AbstractEdgeClassifyConfig; fov_um::NTuple{4,<:Real})
    length(x_um) == length(y_um) ||
        throw(ArgumentError("x_um and y_um must have equal length"))
    fov = (Float64(fov_um[1]), Float64(fov_um[2]), Float64(fov_um[3]), Float64(fov_um[4]))
    fov[1] < fov[2] ||
        throw(ArgumentError("fov_um requires xmin < xmax (got $(fov[1]) >= $(fov[2]))"))
    fov[3] < fov[4] ||
        throw(ArgumentError("fov_um requires ymin < ymax (got $(fov[3]) >= $(fov[4]))"))
    validate(cfg)
    t0 = time_ns()
    x = collect(Float64, x_um); y = collect(Float64, y_um)
    raw = _classify(x, y, fov, cfg)
    return _build_info(raw, cfg, fov, t0)
end

# ---- the shared multi-cell pipeline ------------------------------------------
#
# cell gate (per-config) → relative-density gate → density-adaptive α →
# un-reflected alpha-shape loops → build_mask (split-to-simple + nesting + debris
# cutoff) → per-emitter labeling against the multi-cell mask. The published mask is
# the alpha-shape of the *observed* cell (no FOV reflection); FOV-truncated
# boundary segments are excluded from the membrane band, so a field-of-view cut is
# never mislabeled as membrane.

# Cell-index selection — the only per-config step.
function _cell_indices(x::Vector{Float64}, y::Vector{Float64},
                         fov::NTuple{4,Float64}, cfg::KdeValleyConfig)
    fp = _kde_valley_footprint(x, y, cfg, fov)
    return findall(fp)
end
function _cell_indices(x::Vector{Float64}, y::Vector{Float64},
                         fov::NTuple{4,Float64}, cfg::OuterPolygonConfig)
    n = length(x)
    X = Matrix{Float64}(undef, 2, n)
    @inbounds for i in 1:n
        X[1, i] = x[i]; X[2, i] = y[i]
    end
    cell_mask = _cell_mask(X, cfg.k_list, cfg.rho_k_thresh)
    return findall(cell_mask)
end

function _classify(x::Vector{Float64}, y::Vector{Float64}, fov::NTuple{4,Float64},
                   cfg::AbstractEdgeClassifyConfig)
    n = length(x)
    n == 0 && throw(ErrorException("$(method_name(cfg)): no emitters to classify"))

    idx0 = _cell_indices(x, y, fov, cfg)
    isempty(idx0) && throw(ErrorException(
        "$(method_name(cfg)): density/footprint gate produced no cell points; " *
        "cloud too sparse or parameters too strict"))

    idx = _relative_core_filter(x, y, idx0, cfg.core_radius_nm / 1000, cfg.core_frac)
    length(idx) >= 3 || throw(ErrorException(
        "$(method_name(cfg)): too few cell points after the relative-density gate " *
        "(core_frac=$(cfg.core_frac)); cloud too sparse"))

    # FOV-edge truncation from the GATED TISSUE (not raw emitters): a noise point on the
    # edge must not mark a side truncated for a cell that doesn't reach it. Used below to
    # exclude FOV-cut boundary segments from membrane labeling.
    sides = _truncated_sides(view(x, idx), view(y, idx), fov, cfg.fov_trunc_tol_nm / 1000)

    Xt = Matrix{Float64}(undef, 2, length(idx))
    @inbounds for (k, i) in enumerate(idx)
        Xt[1, k] = x[i]; Xt[2, k] = y[i]
    end
    # Collapse coincident localizations to distinct coordinates: duplicates don't change
    # the boundary, but ≥alpha_knn copies of a point zero its k-NN spacing (→ adaptive
    # α = 0 → valid triangles dropped → the shape collapses) and degrade the Delaunay.
    # Labeling still runs over every original emitter (in_region below).
    size(Xt, 2) > 1 && (Xt = unique(Xt; dims = 2))

    if cfg.alpha_adaptive
        # Multi-scale α: per triangle, min(local carver, conservative envelope).
        # local = alpha_scale × local k-NN spacing (carves dense concavities, bridges
        # sparse gaps); conservative = max(alpha_nm, alpha_scale × cell-median spacing),
        # a per-cell cap whose alpha_nm floor stays loose enough to hold a diffuse
        # background together yet rejects far-reaching low-density-noise protrusions.
        kdist = _knn_distances(Xt, cfg.alpha_knn)
        conservative_um = max(cfg.alpha_nm / 1000, cfg.alpha_scale * median(kdist))
        loops = _local_alpha_shape_loops(Xt, kdist, cfg.alpha_scale, conservative_um)
    else
        loops = _alpha_shape_loops(Xt, cfg.alpha_nm / 1000)
    end
    isempty(loops) && throw(ErrorException(
        "$(method_name(cfg)): no boundary loops found " *
        "(alpha_adaptive=$(cfg.alpha_adaptive)); check inputs or relax the alpha parameters"))

    cells = build_mask(loops; keep_internal = cfg.keep_internal,
                       min_cell_frac = cfg.min_cell_frac, min_hole_frac = cfg.min_hole_frac)
    isempty(cells) && throw(ErrorException(
        "$(method_name(cfg)): no cell survived the min_cell_frac=$(cfg.min_cell_frac) cutoff"))

    class, inside_outer, dist =
        _label_mask(x, y, cells, cfg.membrane_nm / 1000, fov, cfg.fov_trunc_tol_nm / 1000, sides)

    Xorig = Matrix{Float64}(undef, 2, n)
    @inbounds for i in 1:n
        Xorig[1, i] = x[i]; Xorig[2, i] = y[i]
    end
    loop_diags = _compute_loop_diagnostics(loops, x, y, Xorig, fov, _diag_density_thresh(cfg), cells)

    return (; class, inside_outer, dist, outer_polygon = cells[1].outer,
            loops, cells, sides, loop_diags)
end

# ---- per-emitter labeling against the multi-cell mask ------------------------

# Interior iff inside a cell (and not in a kept hole); membrane iff interior and
# within `membrane_um` of a REAL boundary segment (FOV-cut segments excluded);
# else outside. Threads into a Vector{Bool} (one byte per element → race-free),
# then packs to a BitVector.
function _label_mask(x::Vector{Float64}, y::Vector{Float64}, cells::Vector{CellPolygon},
                     membrane_um::Float64, fov::NTuple{4,Float64}, tol_um::Float64,
                     sides::NamedTuple)
    n = length(x)
    inside = Vector{Bool}(undef, n)
    dist = fill(NaN, n)
    Threads.@threads for i in 1:n
        if in_region(x[i], y[i], cells)
            inside[i] = true
            dist[i] = _mask_membrane_dist(x[i], y[i], cells, fov, tol_um, sides)
        else
            inside[i] = false
        end
    end
    class = Vector{Symbol}(undef, n)
    @inbounds for i in 1:n
        if !inside[i]
            class[i] = :outside
        elseif dist[i] < membrane_um
            class[i] = :membrane
        else
            class[i] = :interior
        end
    end
    return class, BitVector(inside), dist
end

# Distance from (qx,qy) to the nearest real (non-FOV-edge) boundary segment across
# every cell's outer ring and its holes. Inf when the only nearby boundary is a
# FOV cut (→ never membrane there).
function _mask_membrane_dist(qx::Float64, qy::Float64, cells::Vector{CellPolygon},
                             fov::NTuple{4,Float64}, tol_um::Float64, sides::NamedTuple)
    best = Inf
    @inbounds for cell in cells
        d = _dist_to_ring_excl_fov(qx, qy, cell.outer, fov, tol_um, sides)
        d < best && (best = d)
        for h in cell.holes
            dh = _dist_to_ring_excl_fov(qx, qy, h, fov, tol_um, sides)
            dh < best && (best = dh)
        end
    end
    return best
end

function _build_info(raw, cfg::AbstractEdgeClassifyConfig, fov::NTuple{4,Float64}, t0::UInt64)
    class = raw.class
    n = length(class)
    n_out = 0; n_mem = 0; n_int = 0
    @inbounds for c in class
        if c === :outside
            n_out += 1
        elseif c === :membrane
            n_mem += 1
        else
            n_int += 1
        end
    end
    return EdgeClassifyInfo(
        n, class, raw.inside_outer, raw.dist,
        raw.outer_polygon, raw.cells, raw.loops, raw.loop_diags,
        cfg, fov, raw.sides, (time_ns() - t0) / 1e9,
        n_out, n_mem, n_int,
    )
end
