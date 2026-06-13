"""
    classify_emitters(smld::BasicSMLD, cfg::AbstractEdgeClassifyConfig) -> (smld, info)
    classify_emitters(x_um, y_um, cfg::AbstractEdgeClassifyConfig; fov_um) -> info

Classify each emitter as `:outside`, `:membrane`, or `:interior`. The concrete
config type selects the strategy by dispatch (`OuterPolygonConfig`, `KdeValleyConfig`).

The SMLD method follows the package `(out, Info)` convention: it returns the smld
(with the primary class mirrored into `smld.metadata["edge_classify_class"]`) and an
[`EdgeClassifyInfo`](@ref). The coordinate method is the computational core and
returns the `info` directly; `fov_um = (xmin, xmax, ymin, ymax)` in µm.
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
    meta["edge_classify_class"] = String.(info.class)
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

# ---- outer_polygon: the pure core (no Info, no IO, no runtime) ---------------
#
# Faithful transcription of the v1 flow: reflect → build Xfull (originals then
# reflected) → multi-K tissue mask → alpha-shape → point-in-polygon + membrane
# band on ORIGINALS. The orchestration order + Xfull column layout are
# parity-load-bearing.
function _classify_polygon(x::Vector{Float64}, y::Vector{Float64},
                           fov::NTuple{4,Float64}, cfg::OuterPolygonConfig)
    n = length(x)
    sides = _truncated_sides(x, y, fov, cfg.fov_trunc_tol_nm / 1000)
    xfull, yfull, n_reflected = _reflect_emitters(
        x, y, fov, sides, cfg.reflect_radius_nm / 1000)
    Xfull = Matrix{Float64}(undef, 2, length(xfull))
    @inbounds for i in eachindex(xfull)
        Xfull[1, i] = xfull[i]; Xfull[2, i] = yfull[i]
    end

    tmask = _tissue_mask(Xfull, cfg.k_list, cfg.rho_k_thresh)
    Xc = Xfull[:, findall(tmask)]
    loops = _alpha_shape_loops(Xc, cfg.alpha_nm / 1000)
    isempty(loops) && throw(ErrorException(
        "no boundary loops found at alpha_nm=$(cfg.alpha_nm); " *
        "check inputs or relax the alpha threshold"))

    poly = loops[1]
    class, inside_outer, dist = _label(x, y, poly, cfg.membrane_nm / 1000)

    Xorig = Matrix{Float64}(undef, 2, n)
    @inbounds for i in 1:n
        Xorig[1, i] = x[i]; Xorig[2, i] = y[i]
    end
    loop_diags = _compute_loop_diagnostics(loops, x, y, Xorig, fov, _diag_density_thresh(cfg))

    return (; class, inside_outer, dist, loops, poly, sides, n_reflected, loop_diags)
end

_classify(x::Vector{Float64}, y::Vector{Float64}, fov::NTuple{4,Float64},
          cfg::OuterPolygonConfig) = _classify_polygon(x, y, fov, cfg)

# ---- kde_valley: own gate + order; reuses the polygon core on the subset -----
function _classify(x::Vector{Float64}, y::Vector{Float64}, fov::NTuple{4,Float64},
                   cfg::KdeValleyConfig)
    n = length(x)
    fp = _kde_valley_footprint(x, y, cfg)
    any(fp) || throw(ErrorException(
        "kde_valley: KDE-valley + footprint gate produced empty tissue " *
        "(sigma_nm=$(cfg.sigma_nm)); cloud too sparse or sigma too small"))
    idx = findall(fp)

    # Polygon geometry on the footprint subset. Empty k_list = internal density
    # gate OFF (no k-NN; the KDE gate already selected tissue). Shared α=600.
    sub = OuterPolygonConfig(alpha_nm = cfg.alpha_nm, membrane_nm = cfg.membrane_nm,
                             reflect_radius_nm = cfg.reflect_radius_nm,
                             fov_trunc_tol_nm = cfg.fov_trunc_tol_nm,
                             k_list = (), rho_k_thresh = 0.0)
    r = _classify_polygon(x[idx], y[idx], fov, sub)

    class = fill(:outside, n)
    inside_outer = falses(n)
    dist = fill(NaN, n)
    @inbounds for (k, i) in enumerate(idx)
        class[i] = r.class[k]
        inside_outer[i] = r.inside_outer[k]
        dist[i] = r.dist[k]
    end

    _enclosure_fill!(class, x, y, cfg)

    return (; class, inside_outer, dist, loops = r.loops, poly = r.poly,
            sides = r.sides, n_reflected = r.n_reflected, loop_diags = r.loop_diags)
end

# ---- shared helpers ----------------------------------------------------------

# Point-in-polygon + membrane band on the ORIGINAL emitters. Threads into a
# Vector{Bool} (one byte per element → race-free; a BitVector would race on
# shared words at @threads chunk boundaries), then packs to a BitVector.
function _label(x::Vector{Float64}, y::Vector{Float64},
                poly::Vector{NTuple{2,Float64}}, membrane_um::Float64)
    n = length(x)
    inside = Vector{Bool}(undef, n)
    dist = fill(NaN, n)
    Threads.@threads for i in 1:n
        if _point_in_polygon(x[i], y[i], poly)
            inside[i] = true
            dist[i] = _dist_to_polygon(x[i], y[i], poly)
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
        raw.poly, raw.loops, raw.loop_diags,
        cfg, fov, raw.sides, raw.n_reflected, (time_ns() - t0) / 1e9,
        n_out, n_mem, n_int,
    )
end
