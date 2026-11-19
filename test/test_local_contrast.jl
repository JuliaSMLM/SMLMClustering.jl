using SMLMClustering
using SMLMData
using Test
using Random
using Statistics
include(joinpath(@__DIR__, "long", "utils", "helpers.jl"))

@testset "LocalContrastFeature" begin

    @testset "config construction + validation" begin
        cfg = LocalContrastFeature()
        @test cfg isa AbstractStatisticsConfig
        @test cfg.density_k == 200
        @test cfg.background_k == 2000
        @test cfg.use_3d === false
        @test cfg.per_dataset === false

        cfg2 = LocalContrastFeature(density_k = 20, background_k = 200,
                                    per_dataset = true)
        @test cfg2.density_k == 20
        @test cfg2.background_k == 200
        @test cfg2.per_dataset === true

        # background_k must be > density_k
        rng = Xoshiro(20260502)
        pts = [(0.5 * randn(rng), 0.5 * randn(rng), 1) for _ in 1:50]
        smld = _make_2d_smld(pts; n_datasets = 1)
        @test_throws ArgumentError cluster_statistics(
            smld, LocalContrastFeature(density_k = 50, background_k = 50))
        @test_throws ArgumentError cluster_statistics(
            smld, LocalContrastFeature(density_k = 100, background_k = 50))
        @test_throws ArgumentError cluster_statistics(
            smld, LocalContrastFeature(density_k = 0, background_k = 10))
    end

    @testset "shape and pass-through semantics" begin
        rng = Xoshiro(20260502)
        pts = [(0.5 * randn(rng), 0.5 * randn(rng), 1) for _ in 1:200]
        smld = _make_2d_smld(pts; n_datasets = 1)

        cfg = LocalContrastFeature(density_k = 10, background_k = 50)
        smld_out, info = cluster_statistics(smld, cfg)

        @test smld_out === smld    # pass-through (V10)
        @test info isa ClusterStatisticsInfo
        @test info.algorithm === :local_contrast
        @test info.statistic_name === :median_local_contrast
        @test info.n_locs_in == 200
        @test info.elapsed_s >= 0

        c = info.extras[:contrast_per_emitter]
        f = info.extras[:log_density_per_emitter]
        @test c isa Vector{Float64}
        @test f isa Vector{Float64}
        @test length(c) == 200
        @test length(f) == 200
    end

    @testset "small group → NaN feature when n ≤ density_k" begin
        rng = Xoshiro(20260502)
        # 8 points, density_k=10: cannot compute → NaN entries.
        pts = [(0.1 * randn(rng), 0.1 * randn(rng), 1) for _ in 1:8]
        smld = _make_2d_smld(pts; n_datasets = 1)
        cfg = LocalContrastFeature(density_k = 10, background_k = 20)
        _, info = cluster_statistics(smld, cfg)
        @test all(isnan, info.extras[:contrast_per_emitter])
        @test all(isnan, info.extras[:log_density_per_emitter])
        @test isnan(info.statistic)
    end

    @testset "per_dataset isolates feature computation" begin
        rng = Xoshiro(20260502)
        pts = Tuple{Float64, Float64, Int}[]
        # Two datasets at the SAME location with different point counts.
        # Pooled: bg=100 spans the union → sees ds1 + ds2 jointly.
        # Per-dataset: bg=100 is per group → only sees own dataset.
        # The densities differ in this overlap regime.
        for _ in 1:200
            push!(pts, (0.5 * rand(rng), 0.5 * rand(rng), 1))
        end
        for _ in 1:200
            push!(pts, (0.5 * rand(rng), 0.5 * rand(rng), 2))
        end
        smld = _make_2d_smld(pts; n_datasets = 2)

        _, info_pooled = cluster_statistics(
            smld, LocalContrastFeature(density_k = 20, background_k = 100,
                                       per_dataset = false))
        _, info_per = cluster_statistics(
            smld, LocalContrastFeature(density_k = 20, background_k = 100,
                                       per_dataset = true))
        # Log-density at fine scale: pooled sees ~2x the density (both datasets
        # in the same area), so pooled log-densities are ≈ +log(2) higher than
        # per-dataset for the same point.
        f_pool = info_pooled.extras[:log_density_per_emitter]
        f_per = info_per.extras[:log_density_per_emitter]
        @test f_pool != f_per
        # Median pooled log-density should be ~log(2) ≈ 0.69 above per-dataset.
        @test median(filter(isfinite, f_pool)) >
              median(filter(isfinite, f_per)) + 0.4
    end

end  # @testset LocalContrastFeature
