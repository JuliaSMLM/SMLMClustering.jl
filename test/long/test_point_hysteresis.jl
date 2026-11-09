using SMLMClustering
using SMLMData
using Test
using Random
using Statistics
include(joinpath(@__DIR__, "utils", "helpers.jl"))

@testset "PointHysteresisConfig backend (thorough)" begin

    @testset "wires with LocalContrastFeature: blob in uniform field is recovered" begin
        rng = Xoshiro(20260502)
        pts = [(5.0 * rand(rng), 5.0 * rand(rng), 1) for _ in 1:1500]
        blob_start = length(pts) + 1
        append!(pts, _blob(rng, 3.0, 3.0, 0.05, 200))
        blob_end = length(pts)
        smld = _make_2d_smld(pts; n_datasets = 1)

        feat_cfg = LocalContrastFeature(density_k = 20, background_k = 400)
        _, info_f = cluster_statistics(smld, feat_cfg)
        c = info_f.extras[:contrast_per_emitter]
        f = info_f.extras[:log_density_per_emitter]
        f_floor = quantile(filter(isfinite, f), 0.35)
        seed = isfinite.(c) .& isfinite.(f) .& (c .> 0.25) .& (f .> f_floor)
        support = isfinite.(c) .& isfinite.(f) .& (c .> -0.05) .& (f .> f_floor)

        _, info = cluster(smld,
            PointHysteresisConfig(graph_k = 20, min_points = 30);
            seed = seed, support = support)
        @test info.n_clusters >= 1
        @test info.cluster_sizes[1] >= 100
    end

end
