#!/usr/bin/env julia
# Edge / membrane classification CLI for SMLMClustering.EdgeClassify.
#
# Loads a BaGoL-format SMLD JLD2 file, runs `classify_emitters` with the chosen
# config, writes stable artifacts under `<out_dir>/<condition>/<cell>/`.
#
# Usage:
#   julia --project=. scripts/edge_classify.jl \
#       --smld <path.jld2> --condition <COND> --cell <CELL> \
#       --out <out_dir> [--params <params.toml>]
#
# params.toml: a `method` key ("outer_polygon" or "kde_valley") plus lowercase
# config fields, e.g.
#   method = "kde_valley"
#   sigma_nm = 150.0
#   alpha_nm = 600.0

using TOML
using Dates
using JLD2
using SMLMClustering

const _USAGE = """
usage: edge_classify.jl --smld PATH --condition COND --cell CELL --out DIR
                        [--params FILE]
"""

function _parse_args(args)
    out = Dict{String,Any}(
        "smld" => nothing, "condition" => nothing, "cell" => nothing,
        "out" => nothing, "params" => nothing)
    i = 1
    while i <= length(args)
        a = args[i]
        if a in ("--smld", "--condition", "--cell", "--out", "--params")
            i + 1 <= length(args) || error("$a requires a value\n$_USAGE")
            out[a[3:end]] = args[i + 1]; i += 2
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

# Build a config from params.toml. `method` selects the type; remaining keys are
# lowercase config fields. Unknown keys raise a clear error from the constructor.
function _build_config(params_path)
    params_path === nothing && return OuterPolygonConfig()
    raw = TOML.parsefile(params_path)
    method = get(raw, "method", "outer_polygon")
    kw = Dict{Symbol,Any}()
    for (k, v) in raw
        k == "method" && continue
        sym = Symbol(k)
        kw[sym] = sym === :k_list ? Tuple(Int.(v)) : v
    end
    if method == "kde_valley"
        return KdeValleyConfig(; kw...)
    elseif method == "outer_polygon"
        return OuterPolygonConfig(; kw...)
    else
        error("unknown method \"$method\"; use \"outer_polygon\" or \"kde_valley\"")
    end
end

function main()
    opts = _parse_args(ARGS)
    cfg = _build_config(opts["params"])

    @info "loading SMLD" smld_path = opts["smld"]
    smld = JLD2.load(opts["smld"], "smld")

    smld_meta = Dict{String,Any}(
        "smld_path"      => abspath(opts["smld"]),
        "smld_size_bytes"=> filesize(opts["smld"]),
        "smld_mtime_utc" => Dates.format(Dates.unix2datetime(mtime(opts["smld"])),
                                         Dates.dateformat"yyyy-mm-ddTHH:MM:SSZ"),
    )

    x = Float64[e.x for e in smld.emitters]
    y = Float64[e.y for e in smld.emitters]
    fov = (Float64(smld.camera.pixel_edges_x[1]), Float64(smld.camera.pixel_edges_x[end]),
           Float64(smld.camera.pixel_edges_y[1]), Float64(smld.camera.pixel_edges_y[end]))

    out_dir = abspath(opts["out"])
    @info "running classify_emitters" method=method_name(cfg) condition=opts["condition"] cell=opts["cell"]
    info = classify_emitters(x, y, cfg; fov_um = fov)

    leaf = joinpath(out_dir, opts["condition"], opts["cell"])
    write_edge_artifacts(leaf, info, x, y;
                         condition = opts["condition"], cell = opts["cell"],
                         smld_input_meta = smld_meta)

    @info "done" runtime_s = round(info.runtime_s; digits=2) leaf
    for f in ("classified.tsv", "polygon_loops.tsv", "loop_diagnostics.csv",
              "params.json", "manifest.json")
        p = joinpath(leaf, f); println("  ", isfile(p) ? "✓" : "✗", "  ", p)
    end
end

main()
