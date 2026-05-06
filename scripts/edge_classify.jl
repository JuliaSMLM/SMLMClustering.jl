#!/usr/bin/env julia
# Edge / membrane classification CLI for SMLMClustering.EdgeClassify.
#
# Loads a BaGoL-format SMLD JLD2 file, runs `classify_emitters`, writes
# stable artifacts under `<out_dir>/<condition>/<cell>/`. See
# docs/src/edge_classify_interface_v1.md for the contract.
#
# Usage:
#   julia --project=. scripts/edge_classify.jl \
#       --smld <path.jld2> --condition <COND> --cell <CELL> \
#       --out <out_dir> [--params <params.toml>] [--renders]

using TOML
using Dates
using JLD2
using SMLMClustering

const _USAGE = """
usage: edge_classify.jl --smld PATH --condition COND --cell CELL --out DIR
                        [--params FILE] [--renders]
"""

function _parse_args(args)
    out = Dict{String,Any}(
        "smld" => nothing, "condition" => nothing, "cell" => nothing,
        "out" => nothing, "params" => nothing, "renders" => false)
    i = 1
    while i <= length(args)
        a = args[i]
        if a in ("--smld", "--condition", "--cell", "--out", "--params")
            i + 1 <= length(args) || error("$a requires a value\n$_USAGE")
            out[a[3:end]] = args[i + 1]; i += 2
        elseif a == "--renders"
            out["renders"] = true; i += 1
        elseif a in ("-h", "--help")
            print(stdout, _USAGE); exit(0)
        else
            error("unknown arg: $a\n$_USAGE")
        end
    end
    for k in ("smld", "condition", "cell", "out")
        out[k] === nothing && error("--$k is required\n$_USAGE")
    end
    return out
end

function _load_params_toml(path::AbstractString)
    raw = TOML.parsefile(path)
    known = Set(["K_LIST","RHO_K_THRESH","ALPHA_NM","REFLECT_RADIUS_NM",
                 "MEMBRANE_NM","FOV_TRUNC_TOL_NM",
                 "METHOD","GRID_PX_NM","GRID_SMOOTH_NM","GRID_MASK_Q",
                 "GRID_MASK_PEAK_FRAC","GRID_OUTER_BUFFER_NM",
                 "CONCAVITY_METRIC_BUFFER_NM",
                 "MASK_CARVE_SIGMA_UM","MASK_CARVE_K_NOISE",
                 "MASK_CARVE_PIXEL_UM","MASK_CARVE_MIN_COMPONENT_FRAC",
                 "MASK_CARVE_FILL_HOLE_MAX_UM2"])
    unknown = setdiff(keys(raw), known)
    isempty(unknown) || error("unknown params keys: $(collect(unknown))")
    kw = Dict{Symbol,Any}()
    haskey(raw, "K_LIST")           && (kw[:K_LIST] = collect(Int, raw["K_LIST"]))
    haskey(raw, "RHO_K_THRESH")     && (kw[:RHO_K_THRESH] = Float64(raw["RHO_K_THRESH"]))
    haskey(raw, "ALPHA_NM")         && (kw[:ALPHA_NM] = Float64(raw["ALPHA_NM"]))
    haskey(raw, "REFLECT_RADIUS_NM")&& (kw[:REFLECT_RADIUS_NM] = Float64(raw["REFLECT_RADIUS_NM"]))
    haskey(raw, "MEMBRANE_NM")      && (kw[:MEMBRANE_NM] = Float64(raw["MEMBRANE_NM"]))
    haskey(raw, "FOV_TRUNC_TOL_NM") && (kw[:FOV_TRUNC_TOL_NM] = Float64(raw["FOV_TRUNC_TOL_NM"]))
    haskey(raw, "METHOD")           && (kw[:METHOD] = String(raw["METHOD"]))
    haskey(raw, "GRID_PX_NM")       && (kw[:GRID_PX_NM] = Float64(raw["GRID_PX_NM"]))
    haskey(raw, "GRID_SMOOTH_NM")   && (kw[:GRID_SMOOTH_NM] = Float64(raw["GRID_SMOOTH_NM"]))
    haskey(raw, "GRID_MASK_Q")      && (kw[:GRID_MASK_Q] = Float64(raw["GRID_MASK_Q"]))
    haskey(raw, "GRID_MASK_PEAK_FRAC") &&
        (kw[:GRID_MASK_PEAK_FRAC] = Float64(raw["GRID_MASK_PEAK_FRAC"]))
    haskey(raw, "GRID_OUTER_BUFFER_NM") &&
        (kw[:GRID_OUTER_BUFFER_NM] = Float64(raw["GRID_OUTER_BUFFER_NM"]))
    haskey(raw, "CONCAVITY_METRIC_BUFFER_NM") &&
        (kw[:CONCAVITY_METRIC_BUFFER_NM] = Float64(raw["CONCAVITY_METRIC_BUFFER_NM"]))
    haskey(raw, "MASK_CARVE_SIGMA_UM") &&
        (kw[:MASK_CARVE_SIGMA_UM] = Float64(raw["MASK_CARVE_SIGMA_UM"]))
    haskey(raw, "MASK_CARVE_K_NOISE") &&
        (kw[:MASK_CARVE_K_NOISE] = Float64(raw["MASK_CARVE_K_NOISE"]))
    haskey(raw, "MASK_CARVE_PIXEL_UM") &&
        (kw[:MASK_CARVE_PIXEL_UM] = Float64(raw["MASK_CARVE_PIXEL_UM"]))
    haskey(raw, "MASK_CARVE_MIN_COMPONENT_FRAC") &&
        (kw[:MASK_CARVE_MIN_COMPONENT_FRAC] = Float64(raw["MASK_CARVE_MIN_COMPONENT_FRAC"]))
    haskey(raw, "MASK_CARVE_FILL_HOLE_MAX_UM2") &&
        (kw[:MASK_CARVE_FILL_HOLE_MAX_UM2] = Float64(raw["MASK_CARVE_FILL_HOLE_MAX_UM2"]))
    return EdgeClassifyParams(; kw...)
end

function main()
    opts = _parse_args(ARGS)
    params = opts["params"] === nothing ? EdgeClassifyParams() :
                                          _load_params_toml(opts["params"])

    @info "loading SMLD" smld_path = opts["smld"]
    smld = JLD2.load(opts["smld"], "smld")

    smld_meta = Dict{String,Any}(
        "smld_path"      => abspath(opts["smld"]),
        "smld_size_bytes"=> filesize(opts["smld"]),
        "smld_mtime_utc" => Dates.format(Dates.unix2datetime(mtime(opts["smld"])),
                                         Dates.dateformat"yyyy-mm-ddTHH:MM:SSZ"),
    )

    out_dir = abspath(opts["out"])
    @info "running classify_emitters" out_dir condition=opts["condition"] cell=opts["cell"]
    result = classify_emitters(smld;
        params = params,
        out_dir = out_dir,
        condition = opts["condition"],
        cell = opts["cell"],
        write_artifacts = true,
        write_renders = opts["renders"],
        smld_input_meta = smld_meta,
    )

    leaf = joinpath(out_dir, opts["condition"], opts["cell"])
    @info "done" runtime_s = round(result.runtime_s; digits=2) leaf
    for f in ("classified.tsv","polygon_loops.tsv","loop_diagnostics.csv",
              "params.json","manifest.json")
        p = joinpath(leaf, f); println("  ", isfile(p) ? "✓" : "✗", "  ", p)
    end
end

main()
