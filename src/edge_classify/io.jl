"""
Internal: artifact writers (TSV/CSV/JSON) and the manifest, driven by an
`EdgeClassifyInfo`. `write_edge_artifacts(leaf, info, x_um, y_um; condition, cell)`
is the public entry; per-config serialization goes through the `to_dict` /
`method_name` traits. JSON is written by hand (no JSON dep).
"""

# ---------- JSON helpers ------------------------------------------------------

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

# ---------- schema versions ---------------------------------------------------

# classified.tsv schema 2: adds the `in_cell` column (topological membership).
const _SCHEMA_VERSION_CLASSIFIED = 2
const _SCHEMA_VERSION_POLYGON_LOOPS = 1
const _SCHEMA_VERSION_LOOP_DIAGNOSTICS = 2
# params.json schema 2: per-config param set (method-specific keys; only the
# fields that ran), produced by `to_dict`.
const _SCHEMA_VERSION_PARAMS = 2
const _SCHEMA_VERSION_MANIFEST = 1
const _SCHEMA_VERSION_EFFECTIVE_OUTER = 1
const _SCHEMA_VERSION_MASK_CARVE_DIAGNOSTIC = 1

# ---------- artifact writers --------------------------------------------------

function _write_classified_tsv(path::AbstractString, info::EdgeClassifyInfo,
                               x_um::Vector{Float64}, y_um::Vector{Float64};
                               condition::AbstractString, cell::AbstractString)
    open(path, "w") do io
        println(io, "# schema_version: ", _SCHEMA_VERSION_CLASSIFIED)
        println(io, "# condition: ", condition)
        println(io, "# cell: ", cell)
        println(io, "# n_emitters: ", info.n_emitters)
        println(io, "# coord_units: um")
        println(io, "emitter_id\tx_um\ty_um\tclass\tinside_outer\tin_cell\tdist_to_outer_um")
        for i in 1:info.n_emitters
            d = info.dist_to_outer_um[i]
            d_str = isnan(d) ? "NaN" : string(round(d; digits=6))
            println(io, i, "\t",
                    round(x_um[i]; digits=6), "\t",
                    round(y_um[i]; digits=6), "\t",
                    String(info.class[i]), "\t",
                    Int(info.inside_outer[i]), "\t",
                    Int(info.class[i] != :outside), "\t",
                    d_str)
        end
    end
end

function _write_polygon_loops_tsv(path::AbstractString, info::EdgeClassifyInfo)
    open(path, "w") do io
        println(io, "# schema_version: ", _SCHEMA_VERSION_POLYGON_LOOPS)
        println(io, "# alpha_nm: ", round(Int, info.config.alpha_nm))
        println(io, "# reflect_radius_nm: ", round(Int, info.config.reflect_radius_nm))
        println(io, "# loop_count: ", length(info.loops))
        println(io, "loop_id\tvertex_id\tx_um\ty_um")
        for (lid, verts) in enumerate(info.loops)
            for (vid, (vx, vy)) in enumerate(verts)
                println(io, lid, "\t", vid, "\t",
                        round(vx; digits=6), "\t",
                        round(vy; digits=6))
            end
        end
    end
end

function _write_loop_diagnostics_csv(path::AbstractString, info::EdgeClassifyInfo)
    open(path, "w") do io
        println(io, "# schema_version: ", _SCHEMA_VERSION_LOOP_DIAGNOSTICS)
        println(io, "loop_id,vertex_count,area_um2,n_emitters_inside,",
                    "frac_in_fov,frac_dense,median_rhoK,used_in_outer,heuristic_type")
        for d in info.loop_diagnostics
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

function _write_effective_outer_tsv(path::AbstractString, info::EdgeClassifyInfo)
    open(path, "w") do io
        println(io, "# schema_version: ", _SCHEMA_VERSION_EFFECTIVE_OUTER)
        println(io, "# method: ", method_name(info.config))
        println(io, "# vertex_count: ", length(info.outer_polygon))
        println(io, "vertex_id\tx_um\ty_um\tmethod")
        for (vid, (vx, vy)) in enumerate(info.outer_polygon)
            println(io, vid, "\t",
                    round(vx; digits=6), "\t",
                    round(vy; digits=6), "\t",
                    method_name(info.config))
        end
    end
end

function _write_mask_carve_diagnostic_json(path::AbstractString, info::EdgeClassifyInfo)
    diag = info.mask_carve_diagnostic
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

function _git_provenance(repo_root::AbstractString)
    sha = ""; clean = true
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

function _write_params_json(path::AbstractString, info::EdgeClassifyInfo;
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
        "n_emitters"       => info.n_emitters,
        "n_reflected"      => info.n_reflected,
        "fov_um"           => collect(Float64, info.fov_um),
        "truncated_sides"  => Dict{String,Any}(
            "L" => info.truncated_sides.L, "R" => info.truncated_sides.R,
            "B" => info.truncated_sides.B, "T" => info.truncated_sides.T),
        "params"           => to_dict(info.config),
        "runtime_s"        => info.runtime_s,
    )
    open(path, "w") do io
        _json_write(io, payload; indent = 0)
        println(io)
    end
end

function _write_manifest_json(path::AbstractString, leaf_dir::AbstractString,
                              out_dir::AbstractString;
                              condition::AbstractString, cell::AbstractString,
                              renders_written::Bool,
                              is_mask_carve::Bool, mask_carve_applied::Bool)
    artifacts = Dict{String,Any}(
        "classified_tsv"       => Dict{String,Any}("path" => "classified.tsv",
                                                    "schema_version" => _SCHEMA_VERSION_CLASSIFIED),
        "polygon_loops_tsv"    => Dict{String,Any}("path" => "polygon_loops.tsv",
                                                    "schema_version" => _SCHEMA_VERSION_POLYGON_LOOPS),
        "loop_diagnostics_csv" => Dict{String,Any}("path" => "loop_diagnostics.csv",
                                                    "schema_version" => _SCHEMA_VERSION_LOOP_DIAGNOSTICS),
        "params_json"          => Dict{String,Any}("path" => "params.json",
                                                    "schema_version" => _SCHEMA_VERSION_PARAMS),
        "classified_png"       => Dict{String,Any}("path" => "classified.png",
                                                    "written" => renders_written,
                                                    "schema_version" => nothing),
        "loop_overlay_png"     => Dict{String,Any}("path" => "loop_overlay.png",
                                                    "written" => renders_written,
                                                    "schema_version" => nothing),
        "effective_outer_tsv"  => Dict{String,Any}(
            "path" => "effective_outer.tsv",
            "written" => is_mask_carve,
            "schema_version" => is_mask_carve ? _SCHEMA_VERSION_EFFECTIVE_OUTER : nothing),
        "mask_carve_diagnostic_json" => Dict{String,Any}(
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

"""
    write_edge_artifacts(leaf, info::EdgeClassifyInfo, x_um, y_um;
                         condition, cell, smld_input_meta=nothing, write_renders=false)

Write the diagnostic artifact set (`classified.tsv`, `polygon_loops.tsv`,
`loop_diagnostics.csv`, `params.json`, `manifest.json`, plus `effective_outer.tsv`
/ `mask_carve_diagnostic.json` for `MaskCarveConfig`) into `leaf`. `x_um`/`y_um`
are the original coordinates (the info does not retain them). Compute and IO are
separate: `classify_emitters` does not write artifacts.
"""
function write_edge_artifacts(leaf::AbstractString, info::EdgeClassifyInfo,
                              x_um::Vector{Float64}, y_um::Vector{Float64};
                              condition::AbstractString, cell::AbstractString,
                              smld_input_meta = nothing, write_renders::Bool = false)
    mkpath(leaf)
    _write_classified_tsv(joinpath(leaf, "classified.tsv"), info, x_um, y_um;
                          condition = condition, cell = cell)
    _write_polygon_loops_tsv(joinpath(leaf, "polygon_loops.tsv"), info)
    _write_loop_diagnostics_csv(joinpath(leaf, "loop_diagnostics.csv"), info)
    _write_params_json(joinpath(leaf, "params.json"), info;
                       smld_input_meta = smld_input_meta)

    is_mask_carve = info.config isa MaskCarveConfig
    mask_carve_applied = false
    if is_mask_carve
        _write_effective_outer_tsv(joinpath(leaf, "effective_outer.tsv"), info)
        if info.mask_carve_diagnostic !== nothing
            _write_mask_carve_diagnostic_json(
                joinpath(leaf, "mask_carve_diagnostic.json"), info)
            mask_carve_applied = info.mask_carve_diagnostic.applied
        end
    end

    out_dir = dirname(dirname(leaf))
    _write_manifest_json(joinpath(leaf, "manifest.json"), leaf, out_dir;
                         condition = condition, cell = cell,
                         renders_written = false,
                         is_mask_carve = is_mask_carve,
                         mask_carve_applied = mask_carve_applied)
    return nothing
end
