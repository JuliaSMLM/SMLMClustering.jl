using SMLMClustering
using SMLMData
using Test
using Random
include(joinpath(@__DIR__, "long", "utils", "helpers.jl"))

@testset "Hierarchical backend" begin

    @testset "config construction" begin
        cfg = HierarchicalConfig(cut_threshold = 200.0)
        @test cfg isa AbstractClusterConfig
        @test cfg.cut_threshold == 200.0
        @test cfg.n_clusters === nothing
        @test cfg.linkage === :ward
        @test cfg.min_points == 5
        @test cfg.use_3d === false
        @test cfg.per_dataset === true
        @test cfg.remove_unclustered === false

        cfg2 = HierarchicalConfig(cut_threshold = 150.0, linkage = :single,
                                   min_points = 3, use_3d = true,
                                   per_dataset = false, remove_unclustered = true)
        @test cfg2.linkage === :single
        @test cfg2.min_points == 3
        @test cfg2.use_3d === true
        @test cfg2.per_dataset === false
        @test cfg2.remove_unclustered === true

        cfg3 = HierarchicalConfig(n_clusters = 3, linkage = :ward)
        @test cfg3.cut_threshold === nothing
        @test cfg3.n_clusters == 3
        @test cfg3.linkage === :ward
    end

    @testset "argument validation" begin
        smld = _make_2d_smld([(0.0, 0.0, 1)])
        # cut_threshold must be positive.
        @test_throws ArgumentError cluster(smld, HierarchicalConfig(cut_threshold = 0.0))
        # n_clusters must be ≥ 1.
        @test_throws ArgumentError cluster(smld, HierarchicalConfig(n_clusters = 0))
        # min_points must be ≥ 1.
        @test_throws ArgumentError cluster(smld, HierarchicalConfig(cut_threshold = 100.0, min_points = 0))
        # Linkage must be one of the four valid symbols.
        @test_throws ArgumentError cluster(smld,
            HierarchicalConfig(cut_threshold = 100.0, linkage = :bogus))
        # Exactly one of cut_threshold / n_clusters must be set — neither.
        @test_throws ArgumentError cluster(smld, HierarchicalConfig())
        # Exactly one of cut_threshold / n_clusters must be set — both.
        @test_throws ArgumentError cluster(smld,
            HierarchicalConfig(cut_threshold = 100.0, n_clusters = 2))
    end

    @testset "use_3d on 2D data errors" begin
        smld = _make_2d_smld([(0.0, 0.0, 1), (1.0, 1.0, 1)])
        cfg = HierarchicalConfig(cut_threshold = 100.0, use_3d = true)
        @test_throws ErrorException cluster(smld, cfg)
    end

end
