using SMLMClustering
using SMLMData
using Test

# Dummy concrete subtype used to exercise the abstract `cluster` fallback.
# Defined at top level because @testset scoping disallows struct definitions
# inside its body.
struct _DummyClusterCfg <: AbstractClusterConfig end

@testset "SMLMClustering.jl" begin

    @testset "types — abstract interface" begin
        @test AbstractClusterConfig <: SMLMData.AbstractSMLMConfig
        @test ClusterInfo <: SMLMData.AbstractSMLMInfo
    end

    @testset "ClusterInfo construction" begin
        info = ClusterInfo(100, 85, 15, 4, [20, 25, 18, 22], :dbscan, 0.123)
        @test info.n_locs_in == 100
        @test info.n_clustered == 85
        @test info.n_noise == 15
        @test info.n_clusters == 4
        @test info.cluster_sizes == [20, 25, 18, 22]
        @test sum(info.cluster_sizes) == info.n_clustered
        @test info.algorithm === :dbscan
        @test info.elapsed_s ≈ 0.123
    end

    @testset "cluster() abstract fallback errors" begin
        # The abstract fallback should refuse for any config subtype without a
        # concrete `cluster` method — exercised via a dummy subtype here.
        cam = IdealCamera(1:8, 1:8, 0.1)
        smld = BasicSMLD(SMLMData.Emitter2DFit{Float64}[], cam, 1, 1, Dict{String,Any}())
        @test_throws ErrorException cluster(smld, _DummyClusterCfg())
    end

    include("test_dbscan.jl")
    include("test_hierarchical.jl")
    include("test_voronoi.jl")

end
