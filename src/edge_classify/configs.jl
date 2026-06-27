"""
Config types for `classify_emitters`.

`AbstractEdgeClassifyConfig` is a sibling of the package's `AbstractClusterConfig`
/ `AbstractStatisticsConfig` (all `<: SMLMData.AbstractSMLMConfig`); each concrete
config is dispatched as a method of `classify_emitters`:

- `OuterPolygonConfig` — reflect → multi-K density gate → alpha-shape outer loop →
  point-in-polygon + membrane band.
- `KdeValleyConfig` — adaptive dSTORM density-valley gate (Gaussian-KDE + valley + footprint
  + enclosure); gates on the original cloud, runs the outer-polygon geometry on the
  footprint subset, then folds enclosure.

Struct fields are lowercase (idiomatic Julia); the UPPERCASE `params.json` keys are
produced only at the serialization boundary by `to_dict`.
"""
abstract type AbstractEdgeClassifyConfig <: SMLMData.AbstractSMLMConfig end

"""
    OuterPolygonConfig(; alpha_nm=300, membrane_nm=100, reflect_radius_nm=1500,
                       fov_trunc_tol_nm=150, k_list=(16,128), rho_k_thresh=200)

Point-in-polygon vs the alpha-shape outer boundary on the FOV-augmented,
multi-K-density-gated set, plus a `membrane_nm` band.
"""
Base.@kwdef struct OuterPolygonConfig <: AbstractEdgeClassifyConfig
    alpha_nm::Float64          = 300.0
    membrane_nm::Float64       = 100.0
    reflect_radius_nm::Float64 = 1500.0
    fov_trunc_tol_nm::Float64  = 150.0
    k_list::Tuple{Vararg{Int}} = (16, 128)   # immutable → provenance-safe
    rho_k_thresh::Float64      = 200.0
end

"""
    KdeValleyConfig(; alpha_nm=600, membrane_nm=100, reflect_radius_nm=1500,
                    fov_trunc_tol_nm=150, sigma_nm=150, ...)

Adaptive density-valley gate for dSTORM data. Gates on the per-FOV KDE density valley
(threshold-free, no per-cell tuning). The defaults are tuned for dSTORM membrane data
— notably `alpha_nm = 600` (vs. the polygon default of 300) — so a bare
`KdeValleyConfig()` is the intended entry point.
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
    c.enclosure_bin_um > 0 || throw(ArgumentError("enclosure_bin_um must be > 0; got $(c.enclosure_bin_um)"))
    (1 <= c.enclosure_min_hits <= 8) || throw(ArgumentError("enclosure_min_hits must be in 1:8; got $(c.enclosure_min_hits)"))
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

# Serialization: UPPERCASE on-disk param dict + the write-only METHOD wire value.
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
