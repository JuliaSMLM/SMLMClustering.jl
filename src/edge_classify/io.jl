"""
Internal: artifact writers (TSV/CSV/JSON) and the manifest file.

Schema versions are emitted in file headers and `manifest.json`. JSON is
written by hand (no JSON dep) for the small fixed schemas in §4d/§4e.
"""

# ---------- helpers -----------------------------------------------------------

# Minimal JSON value writer. Handles String, Bool, Integer, Float64, Vector,
# Dict{String,Any}, NamedTuple, Nothing — sufficient for params.json and
# manifest.json schemas.
function _json_write(io::IO, v; indent::Int=0)
    if v === nothing
        print(io, "null")
    elseif v isa Bool
        print(io, v ? "true" : "false")
    elseif v isa Integer
        print(io, v)
    elseif v isa AbstractFloat
        if isfinite(v)
            print(io, v)
        else
            print(io, "null")
        end
    elseif v isa AbstractString
        _json_write_string(io, v)
    elseif v isa AbstractVector
        print(io, "[")
        for (i, x) in enumerate(v)
            i > 1 && print(io, ", ")
            _json_write(io, x; indent = indent)
        end
        print(io, "]")
    elseif v isa Tuple
        print(io, "[")
        for (i, x) in enumerate(v)
            i > 1 && print(io, ", ")
            _json_write(io, x; indent = indent)
        end
        print(io, "]")
    elseif v isa NamedTuple
        _json_write(io, Dict(string(k) => v[k] for k in keys(v)); indent = indent)
    elseif v isa AbstractDict
        keys_sorted = sort(collect(keys(v)); by = string)
        print(io, "{")
        n = length(keys_sorted)
        for (i, k) in enumerate(keys_sorted)
            print(io, "\n", repeat("  ", indent + 1))
            _json_write_string(io, string(k))
            print(io, ": ")
            _json_write(io, v[k]; indent = indent + 1)
            i < n && print(io, ",")
        end
        n > 0 && print(io, "\n", repeat("  ", indent))
        print(io, "}")
    else
        throw(ArgumentError("unsupported JSON value type: $(typeof(v))"))
    end
end

function _json_write_string(io::IO, s::AbstractString)
    print(io, '"')
    for c in s
        if c == '"'
            print(io, "\\\"")
        elseif c == '\\'
            print(io, "\\\\")
        elseif c == '\n'
            print(io, "\\n")
        elseif c == '\r'
            print(io, "\\r")
        elseif c == '\t'
            print(io, "\\t")
        elseif UInt32(c) < 0x20
            print(io, "\\u", lpad(string(UInt32(c); base=16), 4, '0'))
        else
            print(io, c)
        end
    end
    print(io, '"')
end

# ---------- artifact writers --------------------------------------------------

# schema_version 2 adds the `in_cell` column (topological cell membership;
# differs from `inside_outer` only for METHOD="kde_valley" where enclosure folds
# background into class=="interior" — see EdgeClassificationResult).
const _SCHEMA_VERSION_CLASSIFIED = 2
const _SCHEMA_VERSION_POLYGON_LOOPS = 1
const _SCHEMA_VERSION_LOOP_DIAGNOSTICS = 2
const _SCHEMA_VERSION_PARAMS = 1
# Manifest stays at v1: adding `effective_outer_tsv` and
# `mask_carve_diagnostic_json` artifact entries is additive under the
# existing manifest contract (consumers iterate `artifacts` keys; absence
# / written=false signals "not produced for this method").
const _SCHEMA_VERSION_MANIFEST = 1
const _SCHEMA_VERSION_EFFECTIVE_OUTER = 1
const _SCHEMA_VERSION_MASK_CARVE_DIAGNOSTIC = 1

function _write_classified_tsv(path::AbstractString,
                               result::EdgeClassificationResult,
                               x_um::Vector{Float64}, y_um::Vector{Float64};
                               condition::AbstractString,
                               cell::AbstractString)
    open(path, "w") do io
        println(io, "# schema_version: ", _SCHEMA_VERSION_CLASSIFIED)
        println(io, "# condition: ", condition)
        println(io, "# cell: ", cell)
        println(io, "# n_emitters: ", result.n_emitters)
        println(io, "# coord_units: um")
        println(io, "emitter_id\tx_um\ty_um\tclass\tinside_outer\tin_cell\tdist_to_outer_um")
        for i in 1:result.n_emitters
            d = result.dist_to_outer_um[i]
            d_str = isnan(d) ? "NaN" : string(round(d; digits=6))
            println(io, i, "\t",
                    round(x_um[i]; digits=6), "\t",
                    round(y_um[i]; digits=6), "\t",
                    result.class[i], "\t",
                    Int(result.inside_outer[i]), "\t",
                    Int(result.in_cell[i]), "\t",
                    d_str)
        end
    end
end

function _write_polygon_loops_tsv(path::AbstractString,
                                  result::EdgeClassificationResult)
    open(path, "w") do io
        println(io, "# schema_version: ", _SCHEMA_VERSION_POLYGON_LOOPS)
        println(io, "# alpha_nm: ", round(Int, result.params_used.ALPHA_NM))
        println(io, "# reflect_radius_nm: ", round(Int, result.params_used.REFLECT_RADIUS_NM))
        println(io, "# loop_count: ", length(result.loops))
        println(io, "loop_id\tvertex_id\tx_um\ty_um")
        for (lid, verts) in enumerate(result.loops)
            for (vid, (vx, vy)) in enumerate(verts)
                println(io, lid, "\t", vid, "\t",
                        round(vx; digits=6), "\t",
                        round(vy; digits=6))
            end
        end
    end
end

function _write_loop_diagnostics_csv(path::AbstractString,
                                     result::EdgeClassificationResult)
    open(path, "w") do io
        println(io, "# schema_version: ", _SCHEMA_VERSION_LOOP_DIAGNOSTICS)
        println(io, "loop_id,vertex_count,area_um2,n_emitters_inside,",
                    "frac_in_fov,frac_dense,median_rhoK,used_in_outer,heuristic_type")
        for d in result.loop_diagnostics
            println(io, d.loop_id, ",", d.vertex_count, ",",
                    round(d.area_um2; digits=4), ",",
                    d.n_emitters_inside, ",",
                    round(d.frac_in_fov; digits=4), ",",
                    round(d.frac_dense; digits=4), ",",
                    round(d.median_rhoK; digits=2), ",",
                    d.used_in_outer, ",",
                    d.heuristic_type)
        end
    end
end

function _write_effective_outer_tsv(path::AbstractString,
                                     result::EdgeClassificationResult)
    open(path, "w") do io
        println(io, "# schema_version: ", _SCHEMA_VERSION_EFFECTIVE_OUTER)
        println(io, "# method: ", result.params_used.METHOD)
        println(io, "# vertex_count: ", length(result.outer_polygon))
        println(io, "vertex_id\tx_um\ty_um\tmethod")
        for (vid, (vx, vy)) in enumerate(result.outer_polygon)
            println(io, vid, "\t",
                    round(vx; digits=6), "\t",
                    round(vy; digits=6), "\t",
                    result.params_used.METHOD)
        end
    end
end

function _write_mask_carve_diagnostic_json(path::AbstractString,
                                            result::EdgeClassificationResult)
    diag = result.mask_carve_diagnostic
    payload = Dict{String,Any}(
        "schema_version"             => _SCHEMA_VERSION_MASK_CARVE_DIAGNOSTIC,
        "method"                     => "mask_carve",
        "applied"                    => diag.applied,
        "fallback_reason"            => diag.fallback_reason,
        "params"                     => Dict{String,Any}(
            "MASK_CARVE_SIGMA_UM"           => diag.sigma_um,
            "MASK_CARVE_K_NOISE"            => diag.k_noise,
            "MASK_CARVE_PIXEL_UM"           => diag.pixel_um,
            "MASK_CARVE_MIN_COMPONENT_FRAC" => diag.min_component_frac,
            "MASK_CARVE_FILL_HOLE_MAX_UM2"  => diag.fill_hole_max_um2,
        ),
        "v1_polygon_area_um2"        => diag.v1_polygon_area_um2,
        "carve_polygon_area_um2"     => diag.carve_polygon_area_um2,
        "area_delta_um2"             => diag.area_delta_um2,
        "v1_only_area_um2"           => diag.v1_only_area_um2,
        "carve_only_area_um2"        => diag.carve_only_area_um2,
        "med_v1_carve_distance_um"   => diag.med_v1_carve_distance_um,
        "p95_v1_carve_distance_um"   => diag.p95_v1_carve_distance_um,
        "n_holes_filled"             => diag.n_holes_filled,
        "n_holes_preserved"          => diag.n_holes_preserved,
        "n_carve_polygon_pts"        => diag.n_carve_polygon_pts,
    )
    open(path, "w") do io
        _json_write(io, payload; indent = 0)
        println(io)
    end
end

function _params_to_dict(p::EdgeClassifyConfig)
    return Dict{String,Any}(
        "K_LIST"                        => collect(Int, p.K_LIST),
        "RHO_K_THRESH"                  => p.RHO_K_THRESH,
        "ALPHA_NM"                      => p.ALPHA_NM,
        "REFLECT_RADIUS_NM"             => p.REFLECT_RADIUS_NM,
        "MEMBRANE_NM"                   => p.MEMBRANE_NM,
        "FOV_TRUNC_TOL_NM"              => p.FOV_TRUNC_TOL_NM,
        "METHOD"                        => p.METHOD,
        "GRID_PX_NM"                    => p.GRID_PX_NM,
        "GRID_SMOOTH_NM"                => p.GRID_SMOOTH_NM,
        "GRID_MASK_Q"                   => p.GRID_MASK_Q,
        "GRID_MASK_PEAK_FRAC"           => p.GRID_MASK_PEAK_FRAC,
        "GRID_OUTER_BUFFER_NM"          => p.GRID_OUTER_BUFFER_NM,
        "CONCAVITY_METRIC_BUFFER_NM"    => p.CONCAVITY_METRIC_BUFFER_NM,
        "MASK_CARVE_SIGMA_UM"           => p.MASK_CARVE_SIGMA_UM,
        "MASK_CARVE_K_NOISE"            => p.MASK_CARVE_K_NOISE,
        "MASK_CARVE_PIXEL_UM"           => p.MASK_CARVE_PIXEL_UM,
        "MASK_CARVE_MIN_COMPONENT_FRAC" => p.MASK_CARVE_MIN_COMPONENT_FRAC,
        "MASK_CARVE_FILL_HOLE_MAX_UM2"  => p.MASK_CARVE_FILL_HOLE_MAX_UM2,
        "KDE_SIGMA_NM"                  => p.KDE_SIGMA_NM,
        "KDE_RMAX_SIGMA"                => p.KDE_RMAX_SIGMA,
        "KDE_VALLEY_NBINS"              => p.KDE_VALLEY_NBINS,
        "KDE_VALLEY_FLOORFRAC"          => p.KDE_VALLEY_FLOORFRAC,
        "KDE_VALLEY_SMOOTH"             => p.KDE_VALLEY_SMOOTH,
        "FOOTPRINT_BIN_UM"              => p.FOOTPRINT_BIN_UM,
        "FOOTPRINT_CLOSING_PX"          => p.FOOTPRINT_CLOSING_PX,
        "ENCLOSURE_BIN_UM"              => p.ENCLOSURE_BIN_UM,
        "ENCLOSURE_MIN_HITS"            => p.ENCLOSURE_MIN_HITS,
    )
end

function _git_provenance(repo_root::AbstractString)
    sha = ""
    clean = true
    try
        sha = strip(read(`git -C $repo_root rev-parse HEAD`, String))
    catch
        sha = ""
    end
    try
        s = strip(read(`git -C $repo_root status --porcelain`, String))
        clean = isempty(s)
    catch
        clean = true
    end
    return sha, clean
end

function _write_params_json(path::AbstractString,
                            result::EdgeClassificationResult;
                            smld_input_meta::Union{Nothing,Dict{String,Any}} = nothing,
                            repo_root::AbstractString = pwd())
    sha, clean = _git_provenance(repo_root)
    input_block = smld_input_meta === nothing ?
        Dict{String,Any}("smld_path" => nothing,
                         "smld_mtime_utc" => nothing,
                         "smld_size_bytes" => nothing) :
        smld_input_meta
    payload = Dict{String,Any}(
        "schema_version"   => _SCHEMA_VERSION_PARAMS,
        "git_sha"          => sha,
        "git_status_clean" => clean,
        "timestamp_utc"    => Dates.format(Dates.now(Dates.UTC), Dates.dateformat"yyyy-mm-ddTHH:MM:SSZ"),
        "input"            => input_block,
        "n_emitters"       => result.n_emitters,
        "n_reflected"      => result.n_reflected,
        "fov_um"           => collect(Float64, result.fov_um),
        "truncated_sides"  => Dict{String,Any}(
            "L" => result.truncated_sides.L, "R" => result.truncated_sides.R,
            "B" => result.truncated_sides.B, "T" => result.truncated_sides.T),
        "params"           => _params_to_dict(result.params_used),
        "runtime_s"        => result.runtime_s,
    )
    open(path, "w") do io
        _json_write(io, payload; indent = 0)
        println(io)
    end
end

function _write_manifest_json(path::AbstractString,
                              leaf_dir::AbstractString,
                              out_dir::AbstractString;
                              condition::AbstractString, cell::AbstractString,
                              renders_written::Bool,
                              method::AbstractString = _METHOD_OUTER_POLYGON,
                              mask_carve_applied::Bool = false)
    is_mask_carve = method == _METHOD_MASK_CARVE
    artifacts = Dict{String,Any}(
        "classified_tsv"         => Dict{String,Any}("path" => "classified.tsv",
                                                       "schema_version" => _SCHEMA_VERSION_CLASSIFIED),
        "polygon_loops_tsv"      => Dict{String,Any}("path" => "polygon_loops.tsv",
                                                       "schema_version" => _SCHEMA_VERSION_POLYGON_LOOPS),
        "loop_diagnostics_csv"   => Dict{String,Any}("path" => "loop_diagnostics.csv",
                                                       "schema_version" => _SCHEMA_VERSION_LOOP_DIAGNOSTICS),
        "params_json"            => Dict{String,Any}("path" => "params.json",
                                                       "schema_version" => _SCHEMA_VERSION_PARAMS),
        "classified_png"         => Dict{String,Any}("path" => "classified.png",
                                                       "written" => renders_written,
                                                       "schema_version" => nothing),
        "loop_overlay_png"       => Dict{String,Any}("path" => "loop_overlay.png",
                                                       "written" => renders_written,
                                                       "schema_version" => nothing),
        # mask_carve artifacts: always present in manifest; written/schema_version
        # reflect actual emission (true only when METHOD == "mask_carve").
        "effective_outer_tsv"          => Dict{String,Any}(
            "path" => "effective_outer.tsv",
            "written" => is_mask_carve,
            "schema_version" => is_mask_carve ? _SCHEMA_VERSION_EFFECTIVE_OUTER : nothing),
        "mask_carve_diagnostic_json"   => Dict{String,Any}(
            "path" => "mask_carve_diagnostic.json",
            "written" => is_mask_carve,
            "applied" => mask_carve_applied,
            "schema_version" => is_mask_carve ? _SCHEMA_VERSION_MASK_CARVE_DIAGNOSTIC : nothing),
    )
    payload = Dict{String,Any}(
        "schema_version" => _SCHEMA_VERSION_MANIFEST,
        "condition"      => condition,
        "cell"           => cell,
        "out_dir"        => abspath(out_dir),
        "leaf_dir"       => abspath(leaf_dir),
        "artifacts"      => artifacts,
        "timestamp_utc"  => Dates.format(Dates.now(Dates.UTC), Dates.dateformat"yyyy-mm-ddTHH:MM:SSZ"),
    )
    open(path, "w") do io
        _json_write(io, payload; indent = 0)
        println(io)
    end
end

# ---------- top-level dispatch (called from classify.jl) ----------------------

function _write_artifacts(leaf::AbstractString,
                          result::EdgeClassificationResult;
                          condition::AbstractString, cell::AbstractString,
                          smld_input_meta = nothing,
                          x_um::Union{Nothing,Vector{Float64}} = nothing,
                          y_um::Union{Nothing,Vector{Float64}} = nothing,
                          write_renders::Bool = false)
    # x_um / y_um are passed through from classify_emitters via a thunked
    # wrapper because EdgeClassificationResult does not retain originals.
    # Caller in classify.jl must pre-supply them.
    x_um === nothing && throw(ArgumentError("_write_artifacts requires x_um"))
    y_um === nothing && throw(ArgumentError("_write_artifacts requires y_um"))

    _write_classified_tsv(joinpath(leaf, "classified.tsv"), result, x_um, y_um;
                          condition = condition, cell = cell)
    _write_polygon_loops_tsv(joinpath(leaf, "polygon_loops.tsv"), result)
    _write_loop_diagnostics_csv(joinpath(leaf, "loop_diagnostics.csv"), result)
    _write_params_json(joinpath(leaf, "params.json"), result;
                       smld_input_meta = smld_input_meta)

    is_mask_carve = result.params_used.METHOD == _METHOD_MASK_CARVE
    mask_carve_applied = false
    if is_mask_carve
        _write_effective_outer_tsv(joinpath(leaf, "effective_outer.tsv"), result)
        if result.mask_carve_diagnostic !== nothing
            _write_mask_carve_diagnostic_json(
                joinpath(leaf, "mask_carve_diagnostic.json"), result)
            mask_carve_applied = result.mask_carve_diagnostic.applied
        end
    end

    out_dir = dirname(dirname(leaf))   # leaf == <out>/<cond>/<cell>
    _write_manifest_json(joinpath(leaf, "manifest.json"), leaf, out_dir;
                         condition = condition, cell = cell,
                         renders_written = false,
                         method = result.params_used.METHOD,
                         mask_carve_applied = mask_carve_applied)

    # Renders are deferred — implemented externally for now.
    # write_renders=true currently has no effect inside the package.
    return nothing
end
