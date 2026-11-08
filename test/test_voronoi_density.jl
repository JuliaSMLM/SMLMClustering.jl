using SMLMClustering
using SMLMData
using Test
using Random
using Statistics
include(joinpath(@__DIR__, "long", "utils", "helpers.jl"))

@testset "VoronoiDensity backend" begin

    @testset "config construction" begin
        cfg = VoronoiDensityConfig()
        @test cfg isa AbstractStatisticsConfig
        @test cfg.use_3d === false
        @test cfg.per_dataset === true

        cfg2 = VoronoiDensityConfig(use_3d = false, per_dataset = false)
        @test cfg2.per_dataset === false
    end

    @testset "duplicate (x,y) coordinates raise ArgumentError" begin
        # Exact-coincident generators cause DelaunayTriangulation.get_area
        # to raise KeyError; we want a clean ArgumentError before triangulation.
        pts = [(0.0, 0.0, 1), (1.0, 0.0, 1), (0.0, 1.0, 1), (1.0, 1.0, 1),
               (0.0, 0.0, 1)]  # duplicate of first point
        smld = _make_2d_smld(pts; n_datasets = 1)
        cfg = VoronoiDensityConfig(per_dataset = false)
        @test_throws ArgumentError cluster_statistics(smld, cfg)

        # Per-dataset path also surfaces the error per-group.
        pts2 = [(0.0, 0.0, 1), (1.0, 0.0, 1), (0.0, 1.0, 1), (0.0, 0.0, 1)]
        smld2 = _make_2d_smld(pts2; n_datasets = 1)
        cfg2 = VoronoiDensityConfig(per_dataset = true)
        @test_throws ArgumentError cluster_statistics(smld2, cfg2)
    end

    @testset "use_3d=true raises ArgumentError" begin
        smld = _make_2d_smld([(0.0, 0.0, 1), (1.0, 0.0, 1), (0.0, 1.0, 1)])
        @test_throws ArgumentError cluster_statistics(smld, VoronoiDensityConfig(use_3d = true))
    end

    @testset "passthrough: SMLD reference unchanged, emitter ids untouched" begin
        rng = Xoshiro(7)
        pts = Tuple{Float64,Float64,Int}[]
        for _ in 1:80
            push!(pts, (2.0 * rand(rng), 2.0 * rand(rng), 1))
        end
        smld = _make_2d_smld(pts; n_datasets = 1)
        # Pre-flight: all ids are 0 (default emitter state).
        @test all(e -> e.id == 0, smld.emitters)

        cfg = VoronoiDensityConfig(per_dataset = false)
        smld_out, _ = cluster_statistics(smld, cfg)
        @test smld_out === smld
        @test all(e -> e.id == 0, smld.emitters)
    end

end
