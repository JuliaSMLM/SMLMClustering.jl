"""
Config types for `classify_emitters`.

`AbstractEdgeClassifyConfig` is a sibling of the package's `AbstractClusterConfig`
/ `AbstractStatisticsConfig` (all `<: SMLMData.AbstractSMLMConfig`); each concrete
config is dispatched as a method of `classify_emitters`:

- `OuterPolygonConfig` — multi-K density gate → multi-cell alpha-shape mask →
  point-in-region + membrane band.
- `KdeValleyConfig` — adaptive dSTORM density-valley gate (Gaussian-KDE + valley +
  footprint); gates on the original cloud, then builds the multi-cell mask on the
  footprint subset.

Struct fields are lowercase (idiomatic Julia); the UPPERCASE `params.json` keys are
produced only at the serialization boundary by `to_dict`.
"""
abstract type AbstractEdgeClassifyConfig <: SMLMData.AbstractSMLMConfig end

"""
    OuterPolygonConfig(; alpha_nm=300, membrane_nm=100,
                       fov_trunc_tol_nm=150, k_list=(16,128), rho_k_thresh=200)

Point-in-region vs the multi-cell alpha-shape mask on the multi-K-density-gated set,
plus a `membrane_nm` band.
"""
Base.@kwdef struct OuterPolygonConfig <: AbstractEdgeClassifyConfig
    alpha_nm::Float64          = 300.0
    membrane_nm::Float64       = 100.0
    fov_trunc_tol_nm::Float64  = 150.0
    k_list::Tuple{Vararg{Int}} = (16, 128)   # immutable → provenance-safe
    rho_k_thresh::Float64      = 200.0
    # multi-cell mask pipeline (shared with KdeValleyConfig)
    core_frac::Float64         = 0.10        # relative-density gate (0 disables)
    core_radius_nm::Float64    = 600.0       # gate neighborhood radius
    alpha_adaptive::Bool       = true        # multi-scale α (local carver ∩ per-cell envelope)
    alpha_knn::Int             = 5           # k for the adaptive-α length scale
    alpha_scale::Float64       = 2.0         # ×local k-NN (carver) and ×cell-median (envelope)
    keep_internal::Bool        = false       # keep internal holes (else solid cells)
    min_cell_frac::Float64     = 1/3         # drop cells < frac × largest (0 keeps all)
    min_hole_frac::Float64     = 0.0         # drop holes < frac × cell outer area (0 keeps all)
end

"""
    KdeValleyConfig(; alpha_nm=600, membrane_nm=100,
                    fov_trunc_tol_nm=150, sigma_nm=150, ...)

Adaptive density-valley gate for dSTORM data. Gates on the per-FOV KDE density valley
(threshold-free, no per-cell tuning). The defaults are tuned for dSTORM membrane data
— notably `alpha_nm = 600` (vs. the polygon default of 300) — so a bare
`KdeValleyConfig()` is the intended entry point.
"""
Base.@kwdef struct KdeValleyConfig <: AbstractEdgeClassifyConfig
    alpha_nm::Float64          = 600.0
    membrane_nm::Float64       = 100.0
    fov_trunc_tol_nm::Float64  = 150.0
    sigma_nm::Float64          = 150.0
    rmax_sigma::Float64        = 3.0
    valley_nbins::Int          = 140
    valley_floorfrac::Float64  = 0.05
    valley_smooth::Int         = 4
    footprint_bin_um::Float64  = 0.2
    footprint_closing_px::Int  = 3
    # multi-cell mask pipeline (shared with OuterPolygonConfig)
    core_frac::Float64         = 0.10        # relative-density gate (0 disables)
    core_radius_nm::Float64    = 600.0       # gate neighborhood radius
    alpha_adaptive::Bool       = true        # multi-scale α (local carver ∩ per-cell envelope)
    alpha_knn::Int             = 5           # k for the adaptive-α length scale
    alpha_scale::Float64       = 2.0         # ×local k-NN (carver) and ×cell-median (envelope)
    keep_internal::Bool        = false       # keep internal holes (else solid cells)
    min_cell_frac::Float64     = 1/3         # drop cells < frac × largest (0 keeps all)
    min_hole_frac::Float64     = 0.0         # drop holes < frac × cell outer area (0 keeps all)
end

# ---- validation (per type; called once at dispatch entry) --------------------

function _validate_geom(c::AbstractEdgeClassifyConfig)
    c.alpha_nm > 0           || throw(ArgumentError("alpha_nm must be > 0; got $(c.alpha_nm)"))
    c.membrane_nm >= 0       || throw(ArgumentError("membrane_nm must be >= 0; got $(c.membrane_nm)"))
    c.fov_trunc_tol_nm >= 0  || throw(ArgumentError("fov_trunc_tol_nm must be >= 0; got $(c.fov_trunc_tol_nm)"))
    (0 <= c.core_frac <= 1)  || throw(ArgumentError("core_frac must be in [0,1]; got $(c.core_frac)"))
    c.core_radius_nm > 0     || throw(ArgumentError("core_radius_nm must be > 0; got $(c.core_radius_nm)"))
    c.alpha_knn >= 1         || throw(ArgumentError("alpha_knn must be >= 1; got $(c.alpha_knn)"))
    c.alpha_scale > 0        || throw(ArgumentError("alpha_scale must be > 0; got $(c.alpha_scale)"))
    (0 <= c.min_cell_frac < 1) || throw(ArgumentError("min_cell_frac must be in [0,1); got $(c.min_cell_frac)"))
    (0 <= c.min_hole_frac < 1) || throw(ArgumentError("min_hole_frac must be in [0,1); got $(c.min_hole_frac)"))
    return nothing
end

# Fallback: an unsupported config subtype gets the clear "use a concrete config"
# error here (validate runs before dispatch), mirroring cluster()'s fallback.
validate(c::AbstractEdgeClassifyConfig) =
    error("classify_emitters has no method for config type $(typeof(c)); " *
          "use a concrete config: OuterPolygonConfig or KdeValleyConfig.")

function validate(c::OuterPolygonConfig)
    _validate_geom(c)
    isempty(c.k_list) && throw(ArgumentError("k_list must be non-empty"))
    all(>(0), c.k_list) || throw(ArgumentError("k_list entries must be > 0; got $(c.k_list)"))
    c.rho_k_thresh >= 0 || throw(ArgumentError("rho_k_thresh must be >= 0; got $(c.rho_k_thresh)"))
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
    return nothing
end

# ---- traits ------------------------------------------------------------------

"""
    method_name(cfg) -> String

Short identifier for the edge-classification strategy selected by `cfg`'s concrete
type: `"outer_polygon"` for [`OuterPolygonConfig`](@ref) and `"kde_valley"` for
[`KdeValleyConfig`](@ref). Used in diagnostics, logging, and output manifests.
"""
method_name(::OuterPolygonConfig) = "outer_polygon"
method_name(::KdeValleyConfig)    = "kde_valley"

# Density threshold for the loop-diagnostic `frac_dense` column (OuterPolygonConfig
# only; kde runs the polygon core on a footprint subset via an OuterPolygonConfig).
_diag_density_thresh(c::OuterPolygonConfig) = c.rho_k_thresh
_diag_density_thresh(::KdeValleyConfig) = 0.0   # KDE gate selects cell; no rho_k

# Serialization: UPPERCASE on-disk param dict + the write-only METHOD wire value.
function _geom_dict(c::AbstractEdgeClassifyConfig)
    return Dict{String,Any}(
        "METHOD"            => method_name(c),
        "ALPHA_NM"          => c.alpha_nm,
        "MEMBRANE_NM"       => c.membrane_nm,
        "FOV_TRUNC_TOL_NM"  => c.fov_trunc_tol_nm,
        "CORE_FRAC"         => c.core_frac,
        "CORE_RADIUS_NM"    => c.core_radius_nm,
        "ALPHA_ADAPTIVE"    => c.alpha_adaptive,
        "ALPHA_KNN"         => c.alpha_knn,
        "ALPHA_SCALE"       => c.alpha_scale,
        "KEEP_INTERNAL"     => c.keep_internal,
        "MIN_CELL_FRAC"     => c.min_cell_frac,
        "MIN_HOLE_FRAC"     => c.min_hole_frac,
    )
end

function to_dict(c::OuterPolygonConfig)
    d = _geom_dict(c)
    d["K_LIST"] = collect(Int, c.k_list)
    d["RHO_K_THRESH"] = c.rho_k_thresh
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
    return d
end
