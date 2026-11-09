using SMLMClustering
using SMLMData
using Test
using Random
include(joinpath(@__DIR__, "long", "utils", "helpers.jl"))

@testset "Voronoi backend" begin

    @testset "config construction" begin
        cfg = VoronoiConfig()
        @test cfg isa AbstractClusterConfig
        @test cfg.density_factor == 2.0
        @test cfg.min_points == 5
        @test cfg.use_3d === false
        @test cfg.per_dataset === true
        @test cfg.remove_unclustered === false

        cfg2 = VoronoiConfig(density_factor = 3.5, min_points = 3,
                              per_dataset = false, remove_unclustered = true)
        @test cfg2.density_factor == 3.5
        @test cfg2.min_points == 3
        @test cfg2.per_dataset === false
        @test cfg2.remove_unclustered === true
    end

    @testset "argument validation" begin
        smld = _make_2d_smld([(0.0, 0.0, 1), (1.0, 0.0, 1), (0.0, 1.0, 1)])
        @test_throws ArgumentError cluster(smld, VoronoiConfig(density_factor = 0.0))
        @test_throws ArgumentError cluster(smld, VoronoiConfig(density_factor = -1.0))
        @test_throws ArgumentError cluster(smld, VoronoiConfig(min_points = 0))
        @test_throws ArgumentError cluster(smld, VoronoiConfig(use_3d = true))
    end

    @testset "duplicate coordinates raise ArgumentError" begin
        # Exact-coincident (x,y) generators cause DelaunayTriangulation.get_area
        # to raise KeyError; we want a clean ArgumentError before triangulation.
        pts = [(0.0, 0.0, 1), (1.0, 0.0, 1), (0.0, 1.0, 1), (1.0, 1.0, 1),
               (0.0, 0.0, 1)]  # duplicate of first point
        smld = _make_2d_smld(pts; n_datasets = 1)
        cfg = VoronoiConfig(density_factor = 2.0, min_points = 1, per_dataset = false)
        @test_throws ArgumentError cluster(smld, cfg)

        # Duplicate in a multi-dataset split: per_dataset=true means each dataset
        # is processed independently, so the error fires per-group.
        pts2 = [(0.0, 0.0, 1), (1.0, 0.0, 1), (0.0, 1.0, 1), (0.0, 0.0, 1)]
        smld2 = _make_2d_smld(pts2; n_datasets = 1)
        cfg2 = VoronoiConfig(density_factor = 2.0, min_points = 1, per_dataset = true)
        @test_throws ArgumentError cluster(smld2, cfg2)
    end

end
