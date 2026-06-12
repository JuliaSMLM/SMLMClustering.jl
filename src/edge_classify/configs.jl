"""
Config types for `classify_emitters`.

`AbstractEdgeClassifyConfig` is a sibling of the package's `AbstractClusterConfig`
/ `AbstractStatisticsConfig` (all `<: SMLMData.AbstractSMLMConfig`); each concrete
config is dispatched as a method of `classify_emitters`. `AbstractPolygonConfig`
is the shared family — reflect → multi-K density gate → alpha-shape outer loop →
point-in-polygon + membrane band — whose members differ only in two dispatched
hooks (`_effective_polygon`, `_refine!`). `KdeValleyConfig` is its own family
(gates on the original cloud, runs the polygon core on a footprint subset, folds
enclosure) and carries the validated `alpha_nm = 600` as its own default.

Struct fields are lowercase (idiomatic Julia); the UPPERCASE `params.json` keys
are produced only at the serialization boundary by `to_dict`.
"""

abstract type AbstractEdgeClassifyConfig <: SMLMData.AbstractSMLMConfig end

abstract type AbstractPolygonConfig <: AbstractEdgeClassifyConfig end

"""
    OuterPolygonConfig(; alpha_nm=300, membrane_nm=100, reflect_radius_nm=1500,
                       fov_trunc_tol_nm=150, k_list=[16,128], rho_k_thresh=200)

v1 outer-polygon classifier: point-in-polygon vs the alpha-shape outer boundary
on the FOV-augmented, multi-K-density-gated set, plus a `membrane_nm` band.
"""
Base.@kwdef struct OuterPolygonConfig <: AbstractPolygonConfig
    alpha_nm::Float64          = 300.0
    membrane_nm::Float64       = 100.0
    reflect_radius_nm::Float64 = 1500.0
    fov_trunc_tol_nm::Float64  = 150.0
    k_list::Vector{Int}        = [16, 128]
    rho_k_thresh::Float64      = 200.0
end

"""
    GridHybridConfig(; <OuterPolygon fields>, grid_px_nm=50, grid_smooth_nm=80,
                     grid_mask_q=0.03, grid_mask_peak_frac=0.26, grid_outer_buffer_nm=800)

Outer-polygon topology + a density-grid post-pass that promotes interior emitters
near the boundary to `membrane`. Never changes `outside`.
"""
Base.@kwdef struct GridHybridConfig <: AbstractPolygonConfig
    alpha_nm::Float64             = 300.0
    membrane_nm::Float64          = 100.0
    reflect_radius_nm::Float64    = 1500.0
    fov_trunc_tol_nm::Float64     = 150.0
    k_list::Vector{Int}           = [16, 128]
    rho_k_thresh::Float64         = 200.0
    grid_px_nm::Float64           = 50.0
    grid_smooth_nm::Float64       = 80.0
    grid_mask_q::Float64          = 0.03
    grid_mask_peak_frac::Float64  = 0.26
    grid_outer_buffer_nm::Float64 = 800.0
end

"""
    MaskCarveConfig(; <OuterPolygon fields>, sigma_um=0.08, k_noise=3.0,
                    pixel_um=0.04, min_component_frac=0.05, fill_hole_max_um2=0.5)

Outer-polygon, but the effective classification polygon is a density-mask carve of
v1 (carve ⊆ v1; never expands outward).
"""
Base.@kwdef struct MaskCarveConfig <: AbstractPolygonConfig
    alpha_nm::Float64           = 300.0
    membrane_nm::Float64        = 100.0
    reflect_radius_nm::Float64  = 1500.0
    fov_trunc_tol_nm::Float64   = 150.0
    k_list::Vector{Int}         = [16, 128]
    rho_k_thresh::Float64       = 200.0
    sigma_um::Float64           = 0.080
    k_noise::Float64            = 3.0
    pixel_um::Float64           = 0.040
    min_component_frac::Float64 = 0.05
    fill_hole_max_um2::Float64  = 0.5
end

"""
    KdeValleyConfig(; alpha_nm=600, membrane_nm=100, reflect_radius_nm=1500,
                    fov_trunc_tol_nm=150, sigma_nm=150, ...)

Validated adaptive dSTORM gate (genmab): continuous Gaussian-KDE density +
background/cell valley threshold + footprint fill, then the outer-polygon geometry
on the footprint subset, then ray-cast enclosure reclass folding enclosed
background into `interior`. Per-FOV adaptive — no per-cell density tuning.
Defaults reproduce the A431 dSTORM validated set; `alpha_nm = 600` is the
validated value (not the polygon-family default 300).
"""
Base.@kwdef struct KdeValleyConfig <: AbstractEdgeClassifyConfig
    alpha_nm::Float64          = 600.0
    membrane_nm::Float64       = 100.0
    reflect_radius_nm::Float64 = 1500.0
    fov_trunc_tol_nm::Float64  = 150.0
    sigma_nm::Float64          = 150.0
    rmax_sigma::Float64        = 3.0
    valley_nbins::Int          = 140
    valley_floorfrac::Float64  = 0.05
    valley_smooth::Int         = 4
    footprint_bin_um::Float64  = 0.2
    footprint_closing_px::Int  = 3
    enclosure_bin_um::Float64  = 0.2
    enclosure_min_hits::Int    = 6
end

# ---- validation (per type; called once at dispatch entry) --------------------

function _validate_geom(c::AbstractEdgeClassifyConfig)
    c.alpha_nm > 0           || throw(ArgumentError("alpha_nm must be > 0; got $(c.alpha_nm)"))
    c.membrane_nm >= 0       || throw(ArgumentError("membrane_nm must be >= 0; got $(c.membrane_nm)"))
    c.reflect_radius_nm >= 0 || throw(ArgumentError("reflect_radius_nm must be >= 0; got $(c.reflect_radius_nm)"))
    c.fov_trunc_tol_nm >= 0  || throw(ArgumentError("fov_trunc_tol_nm must be >= 0; got $(c.fov_trunc_tol_nm)"))
    return nothing
end

function _validate_polygon_gate(c::AbstractPolygonConfig)
    isempty(c.k_list) && throw(ArgumentError("k_list must be non-empty"))
    all(>(0), c.k_list) || throw(ArgumentError("k_list entries must be > 0; got $(c.k_list)"))
    c.rho_k_thresh >= 0 || throw(ArgumentError("rho_k_thresh must be >= 0; got $(c.rho_k_thresh)"))
    return nothing
end

function validate(c::OuterPolygonConfig)
    _validate_geom(c); _validate_polygon_gate(c); return nothing
end

function validate(c::GridHybridConfig)
    _validate_geom(c); _validate_polygon_gate(c)
    c.grid_px_nm > 0     || throw(ArgumentError("grid_px_nm must be > 0; got $(c.grid_px_nm)"))
    c.grid_smooth_nm > 0 || throw(ArgumentError("grid_smooth_nm must be > 0; got $(c.grid_smooth_nm)"))
    (0 <= c.grid_mask_q <= 1)         || throw(ArgumentError("grid_mask_q must be in [0,1]; got $(c.grid_mask_q)"))
    (0 <= c.grid_mask_peak_frac <= 1) || throw(ArgumentError("grid_mask_peak_frac must be in [0,1]; got $(c.grid_mask_peak_frac)"))
    c.grid_outer_buffer_nm >= 0       || throw(ArgumentError("grid_outer_buffer_nm must be >= 0; got $(c.grid_outer_buffer_nm)"))
    return nothing
end

function validate(c::MaskCarveConfig)
    _validate_geom(c); _validate_polygon_gate(c)
    c.sigma_um > 0            || throw(ArgumentError("sigma_um must be > 0; got $(c.sigma_um)"))
    c.k_noise > 0             || throw(ArgumentError("k_noise must be > 0; got $(c.k_noise)"))
    c.pixel_um > 0            || throw(ArgumentError("pixel_um must be > 0; got $(c.pixel_um)"))
    c.min_component_frac >= 0 || throw(ArgumentError("min_component_frac must be >= 0; got $(c.min_component_frac)"))
    c.fill_hole_max_um2 >= 0  || throw(ArgumentError("fill_hole_max_um2 must be >= 0; got $(c.fill_hole_max_um2)"))
    return nothing
end

function validate(c::KdeValleyConfig)
    _validate_geom(c)
    c.sigma_nm > 0   || throw(ArgumentError("sigma_nm must be > 0; got $(c.sigma_nm)"))
    c.rmax_sigma > 0 || throw(ArgumentError("rmax_sigma must be > 0; got $(c.rmax_sigma)"))
    c.valley_nbins > 1 || throw(ArgumentError("valley_nbins must be > 1; got $(c.valley_nbins)"))
    (0 <= c.valley_floorfrac <= 1) || throw(ArgumentError("valley_floorfrac must be in [0,1]; got $(c.valley_floorfrac)"))
    c.valley_smooth >= 0 || throw(ArgumentError("valley_smooth must be >= 0; got $(c.valley_smooth)"))
    c.footprint_bin_um > 0 || throw(ArgumentError("footprint_bin_um must be > 0; got $(c.footprint_bin_um)"))
    c.footprint_closing_px >= 0 || throw(ArgumentError("footprint_closing_px must be >= 0; got $(c.footprint_closing_px)"))
    c.enclosure_bin_um > 0 || throw(ArgumentError("enclosure_bin_um must be > 0; got $(c.enclosure_bin_um)"))
    (1 <= c.enclosure_min_hits <= 8) || throw(ArgumentError("enclosure_min_hits must be in 1:8; got $(c.enclosure_min_hits)"))
    return nothing
end

# ---- traits ------------------------------------------------------------------

method_name(::OuterPolygonConfig) = "outer_polygon"
method_name(::GridHybridConfig)   = "grid_hybrid"
method_name(::MaskCarveConfig)     = "mask_carve"
method_name(::KdeValleyConfig)     = "kde_valley"

# Density threshold for the loop-diagnostic `frac_dense` column. The polygon
# family uses its tissue gate; kde_valley gates by KDE valley (not rho_k), so
# frac_dense is reported against 0 (fraction of vertices with positive density).
_diag_density_thresh(c::AbstractPolygonConfig) = c.rho_k_thresh
_diag_density_thresh(::KdeValleyConfig)        = 0.0

# Serialization: UPPERCASE on-disk param dict + the write-only METHOD wire value.
# Each config records only the fields that actually ran (fixes the v1 wart where
# a kde_valley run logged inert rho_k_thresh/k_list).
function _geom_dict(c::AbstractEdgeClassifyConfig)
    return Dict{String,Any}(
        "METHOD"            => method_name(c),
        "ALPHA_NM"          => c.alpha_nm,
        "MEMBRANE_NM"       => c.membrane_nm,
        "REFLECT_RADIUS_NM" => c.reflect_radius_nm,
        "FOV_TRUNC_TOL_NM"  => c.fov_trunc_tol_nm,
    )
end

function to_dict(c::OuterPolygonConfig)
    d = _geom_dict(c)
    d["K_LIST"] = collect(Int, c.k_list)
    d["RHO_K_THRESH"] = c.rho_k_thresh
    return d
end

function to_dict(c::GridHybridConfig)
    d = _geom_dict(c)
    d["K_LIST"] = collect(Int, c.k_list)
    d["RHO_K_THRESH"] = c.rho_k_thresh
    d["GRID_PX_NM"] = c.grid_px_nm
    d["GRID_SMOOTH_NM"] = c.grid_smooth_nm
    d["GRID_MASK_Q"] = c.grid_mask_q
    d["GRID_MASK_PEAK_FRAC"] = c.grid_mask_peak_frac
    d["GRID_OUTER_BUFFER_NM"] = c.grid_outer_buffer_nm
    return d
end

function to_dict(c::MaskCarveConfig)
    d = _geom_dict(c)
    d["K_LIST"] = collect(Int, c.k_list)
    d["RHO_K_THRESH"] = c.rho_k_thresh
    d["MASK_CARVE_SIGMA_UM"] = c.sigma_um
    d["MASK_CARVE_K_NOISE"] = c.k_noise
    d["MASK_CARVE_PIXEL_UM"] = c.pixel_um
    d["MASK_CARVE_MIN_COMPONENT_FRAC"] = c.min_component_frac
    d["MASK_CARVE_FILL_HOLE_MAX_UM2"] = c.fill_hole_max_um2
    return d
end

function to_dict(c::KdeValleyConfig)
    d = _geom_dict(c)
    d["KDE_SIGMA_NM"] = c.sigma_nm
    d["KDE_RMAX_SIGMA"] = c.rmax_sigma
    d["KDE_VALLEY_NBINS"] = c.valley_nbins
    d["KDE_VALLEY_FLOORFRAC"] = c.valley_floorfrac
    d["KDE_VALLEY_SMOOTH"] = c.valley_smooth
    d["FOOTPRINT_BIN_UM"] = c.footprint_bin_um
    d["FOOTPRINT_CLOSING_PX"] = c.footprint_closing_px
    d["ENCLOSURE_BIN_UM"] = c.enclosure_bin_um
    d["ENCLOSURE_MIN_HITS"] = c.enclosure_min_hits
    return d
end
