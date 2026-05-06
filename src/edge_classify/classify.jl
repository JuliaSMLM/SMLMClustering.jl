"""
    classify_emitters(x_um, y_um; fov_um, params=EdgeClassifyParams(),
                      out_dir=nothing, condition=nothing, cell=nothing,
                      write_artifacts=false, write_renders=false)
        -> EdgeClassificationResult

Classify each (x_um, y_um) emitter as `"outside"`, `"membrane"`, or
`"interior"` using the v1 outer-polygon decision (point-in-polygon vs the
alpha-shape outer boundary on the FOV-augmented multi-K-density-gated set,
plus a `MEMBRANE_NM` band around the boundary).

# Required arguments
- `x_um::AbstractVector{<:Real}`, `y_um::AbstractVector{<:Real}`: original
  emitter coordinates in µm. Same length.
- `fov_um::NTuple{4,Float64}`: camera FOV bounds in µm, ordered
  `(xmin_um, xmax_um, ymin_um, ymax_um)`. Validated `xmin < xmax` and
  `ymin < ymax`.

# Keyword arguments
- `params::EdgeClassifyParams`: pipeline parameters.
- `out_dir, condition, cell`: required when `write_artifacts=true`. Output
  goes to `<out_dir>/<condition>/<cell>/`.
- `write_artifacts::Bool`: emit `classified.tsv`, `polygon_loops.tsv`,
  `loop_diagnostics.csv`, `params.json`, `manifest.json`.
- `write_renders::Bool`: also emit diagnostic PNG renders (currently a
  no-op placeholder; renders live outside the package proper for now).

# Returns
An `EdgeClassificationResult`. Class labels partition the input set.
"""
function classify_emitters(
    x_um::AbstractVector{<:Real},
    y_um::AbstractVector{<:Real};
    fov_um::NTuple{4,Float64},
    params::EdgeClassifyParams = EdgeClassifyParams(),
    out_dir::Union{Nothing,AbstractString} = nothing,
    condition::Union{Nothing,AbstractString} = nothing,
    cell::Union{Nothing,AbstractString} = nothing,
    write_artifacts::Bool = false,
    write_renders::Bool = false,
    smld_input_meta::Union{Nothing,Dict{String,Any}} = nothing,
)
    length(x_um) == length(y_um) ||
        throw(ArgumentError("x_um and y_um must have equal length"))
    fov_um[1] < fov_um[2] ||
        throw(ArgumentError("fov_um requires xmin < xmax (got $(fov_um[1]) >= $(fov_um[2]))"))
    fov_um[3] < fov_um[4] ||
        throw(ArgumentError("fov_um requires ymin < ymax (got $(fov_um[3]) >= $(fov_um[4]))"))
    params.METHOD in _VALID_METHODS ||
        throw(ArgumentError("params.METHOD must be one of $(_VALID_METHODS); got \"$(params.METHOD)\""))
    if params.METHOD == _METHOD_CONCAVE_REFINED
        throw(ArgumentError(
            "METHOD=\"concave_refined\" is reserved for the concave-membrane " *
            "branch and is not implemented yet; use METHOD=\"outer_polygon\", " *
            "\"grid_hybrid\", or \"mask_carve\""))
    end
    if params.METHOD == _METHOD_MASK_CARVE
        params.MASK_CARVE_SIGMA_UM > 0 ||
            throw(ArgumentError("MASK_CARVE_SIGMA_UM must be positive; got $(params.MASK_CARVE_SIGMA_UM)"))
        params.MASK_CARVE_PIXEL_UM > 0 ||
            throw(ArgumentError("MASK_CARVE_PIXEL_UM must be positive; got $(params.MASK_CARVE_PIXEL_UM)"))
        params.MASK_CARVE_K_NOISE > 0 ||
            throw(ArgumentError("MASK_CARVE_K_NOISE must be positive; got $(params.MASK_CARVE_K_NOISE)"))
        params.MASK_CARVE_MIN_COMPONENT_FRAC >= 0 ||
            throw(ArgumentError("MASK_CARVE_MIN_COMPONENT_FRAC must be >= 0; got $(params.MASK_CARVE_MIN_COMPONENT_FRAC)"))
        params.MASK_CARVE_FILL_HOLE_MAX_UM2 >= 0 ||
            throw(ArgumentError("MASK_CARVE_FILL_HOLE_MAX_UM2 must be >= 0; got $(params.MASK_CARVE_FILL_HOLE_MAX_UM2)"))
    end
    if write_artifacts
        out_dir === nothing &&
            throw(ArgumentError("write_artifacts=true requires out_dir"))
        condition === nothing &&
            throw(ArgumentError("write_artifacts=true requires condition"))
        cell === nothing &&
            throw(ArgumentError("write_artifacts=true requires cell"))
    end

    t0 = time()
    x = collect(Float64, x_um); y = collect(Float64, y_um)
    n = length(x)

    # FOV truncation detection + reflection.
    sides = _truncated_sides(x, y, fov_um, params.FOV_TRUNC_TOL_NM / 1000)
    xfull, yfull, n_reflected = _reflect_emitters(
        x, y, fov_um, sides, params.REFLECT_RADIUS_NM / 1000)
    Xfull = Matrix{Float64}(undef, 2, length(xfull))
    @inbounds for i in eachindex(xfull)
        Xfull[1, i] = xfull[i]; Xfull[2, i] = yfull[i]
    end

    # Multi-K density tissue mask on augmented set.
    tmask = _tissue_mask(Xfull, params.K_LIST, params.RHO_K_THRESH)

    # Alpha-shape on tissue points.
    tissue_idx = findall(tmask)
    Xc = Xfull[:, tissue_idx]
    loops = _alpha_shape_loops(Xc, params.ALPHA_NM / 1000)
    isempty(loops) && throw(ErrorException(
        "no boundary loops found at ALPHA_NM=$(params.ALPHA_NM); " *
        "check inputs or relax the alpha threshold"))

    v1_outer = loops[1]
    mask_carve_diag::Union{Nothing, MaskCarveDiagnostic} = nothing
    if params.METHOD == _METHOD_MASK_CARVE
        carve_poly, mask_carve_diag = _build_mask_carve(v1_outer, x, y, fov_um, params)
        outer_polygon = carve_poly
    else
        outer_polygon = v1_outer
    end

    # Per-emitter classification on ORIGINALS only — uses the EFFECTIVE
    # outer polygon (carve for mask_carve when applied; v1 outer otherwise).
    inside_outer = falses(n)
    dist_to_outer = fill(NaN, n)
    Threads.@threads for i in 1:n
        if _point_in_polygon(x[i], y[i], outer_polygon)
            inside_outer[i] = true
        end
    end
    Threads.@threads for i in 1:n
        if inside_outer[i]
            dist_to_outer[i] = _dist_to_polygon(x[i], y[i], outer_polygon)
        end
    end

    membrane_um = params.MEMBRANE_NM / 1000
    class = Vector{String}(undef, n)
    @inbounds for i in 1:n
        if !inside_outer[i]
            class[i] = "outside"
        elseif dist_to_outer[i] < membrane_um
            class[i] = "membrane"
        else
            class[i] = "interior"
        end
    end

    if params.METHOD == _METHOD_GRID_HYBRID
        _apply_grid_hybrid!(class, x, y, dist_to_outer, fov_um, params)
    end

    # Per-loop diagnostics — uses originals-only KDTree.
    Xorig = Matrix{Float64}(undef, 2, n)
    @inbounds for i in 1:n
        Xorig[1, i] = x[i]; Xorig[2, i] = y[i]
    end
    loop_diags = _compute_loop_diagnostics(loops, x, y, Xorig, fov_um, params)

    runtime_s = time() - t0
    result = EdgeClassificationResult(
        n, class, inside_outer, dist_to_outer,
        outer_polygon, loops, loop_diags,
        params, fov_um, sides, n_reflected, runtime_s,
        mask_carve_diag,
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
