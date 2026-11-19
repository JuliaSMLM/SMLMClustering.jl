using SMLMClustering
using SMLMData
using Test
using Random
include(joinpath(@__DIR__, "long", "utils", "helpers.jl"))

@testset "Hopkins backend" begin

    @testset "config construction" begin
        cfg = HopkinsConfig()
        @test cfg isa AbstractStatisticsConfig
        @test cfg.n_samples == 20
        @test cfg.random_repeats == 1
        @test cfg.seed === nothing
        @test cfg.use_3d === false
        @test cfg.per_dataset === true
        @test cfg.region === nothing

        cfg2 = HopkinsConfig(n_samples = 50, random_repeats = 5, seed = 42,
                             use_3d = true, per_dataset = false)
        @test cfg2.n_samples == 50
        @test cfg2.random_repeats == 5
        @test cfg2.seed == 42
        @test cfg2.use_3d === true
        @test cfg2.per_dataset === false
    end

    @testset "input SMLD is not modified (passthrough preserves emitter state)" begin
        rng = Xoshiro(20260427)
        pts = Tuple{Float64,Float64,Int}[]
        for _ in 1:200
            push!(pts, (2.5 + 0.01 * randn(rng), 2.5 + 0.01 * randn(rng), 1))
        end
        smld = _make_2d_smld(pts; n_datasets = 1)
        # Pre-flight: all ids are 0 (default emitter state).
        @test all(e -> e.id == 0, smld.emitters)

        cfg = HopkinsConfig(n_samples = 30, random_repeats = 2, seed = 4,
                            per_dataset = false)
        smld_out, _ = cluster_statistics(smld, cfg)
        @test smld_out === smld
        @test all(e -> e.id == 0, smld.emitters)
    end

    @testset "argument validation" begin
        smld = _make_2d_smld([(1.0, 1.0, 1), (2.0, 2.0, 1), (3.0, 3.0, 1)])
        @test_throws ArgumentError cluster_statistics(smld, HopkinsConfig(n_samples = 0))
        @test_throws ArgumentError cluster_statistics(smld, HopkinsConfig(random_repeats = 0))
        tri = [(0.0, 0.0), (4.0, 0.0), (2.0, 4.0)]
        # region (a 2D polygon) is incompatible with 3D
        @test_throws ArgumentError cluster_statistics(smld, HopkinsConfig(use_3d = true, region = tri))
        # :metadata requested but no polygon present
        @test_throws ArgumentError cluster_statistics(smld, HopkinsConfig(region = :metadata, per_dataset = false))
        # malformed polygon (< 3 vertices)
        @test_throws ArgumentError cluster_statistics(smld, HopkinsConfig(region = [(0.0, 0.0), (1.0, 1.0)], per_dataset = false))
    end

end
