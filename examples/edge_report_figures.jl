# Example + manual coverage for the edge-mask report and figure series.
#
# The figures live in `SMLMClusteringFiguresExt`, so run this with CairoMakie AND
# SMLMRender available (both are deps of this examples project):
#
#     julia --project=examples examples/edge_report_figures.jl
#
# Package extensions are not exercised by `Pkg.test`, so this script is the manual
# smoke-test that `compute_edge_report` / `write_edge_report` / `plot_edge_report`
# stay wired correctly.

using SMLMClustering, SMLMData, CairoMakie, SMLMRender, Random

"""
    synthetic_cell(; n=6000, seed=1) -> BasicSMLD

A synthetic dense, membrane-like cell: a disk with a concave bay carved out of the
right side, plus a little sparse noise — an `Emitter2DFit` SMLD with σ ≈ 20 nm so the
precision-based circle render has something to draw.
"""
function synthetic_cell(; n = 6000, seed = 1)
    rng = MersenneTwister(seed)
    xs = Float64[]; ys = Float64[]
    while length(xs) < n
        x = rand(rng) * 20; y = rand(rng) * 20
        r = hypot(x - 10, y - 10)
        inbay = x > 12 && abs(y - 10) < 3            # carve a bay on the right
        if r < 8 && !inbay
            push!(xs, x); push!(ys, y)
        end
    end
    for _ in 1:50                                     # a little sparse noise outside
        push!(xs, rand(rng) * 20); push!(ys, rand(rng) * 20)
    end
    em = [SMLMData.Emitter2DFit{Float64}(xs[i], ys[i], 1000.0, 5.0, 0.02, 0.02, 30.0, 5.0; id = i)
          for i in eachindex(xs)]
    cam = SMLMData.IdealCamera(1:200, 1:200, 0.1)     # 200×200 px @ 0.1 µm → 20 µm FOV
    return SMLMData.BasicSMLD(em, cam, 1, 1)
end

smld = synthetic_cell()
smld_out, info = classify_emitters(smld, KdeValleyConfig(sigma_nm = 200.0))
@info "classified" n = info.n_emitters interior = info.n_interior membrane = info.n_membrane outside = info.n_outside cells = length(info.cells)

report = compute_edge_report(smld_out, info)
outdir = joinpath(@__DIR__, "..", "dev", "output", "edge_report_example")
write_edge_report(report; output_dir = outdir, condition = "synthetic", cell = "disk")
paths = plot_edge_report(report; output_dir = outdir, prefix = "edge")
@info "wrote edge-mask figures" outdir paths
