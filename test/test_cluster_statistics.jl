using SMLMClustering
using SMLMData
using Test

# Dummy concrete subtype used to exercise the abstract `cluster_statistics`
# fallback. Defined at top level because @testset scoping disallows struct
# definitions inside its body.
struct _DummyStatsCfg <: AbstractStatisticsConfig end

@testset "cluster_statistics — abstract interface" begin

    @testset "type hierarchy" begin
        @test AbstractStatisticsConfig <: SMLMData.AbstractSMLMConfig
        @test ClusterStatisticsInfo <: SMLMData.AbstractSMLMInfo
        @test HopkinsConfig <: AbstractStatisticsConfig
        @test VoronoiDensityConfig <: AbstractStatisticsConfig
    end

    @testset "ClusterStatisticsInfo construction" begin
        info = ClusterStatisticsInfo(
            100, 0.83, :hopkins, :hopkins, 0.012,
            Dict{Symbol,Any}(:hopkins_per_dataset => [0.81, 0.85]),
        )
        @test info.n_locs_in == 100
        @test info.statistic ≈ 0.83
        @test info.statistic_name === :hopkins
        @test info.algorithm === :hopkins
        @test info.elapsed_s ≈ 0.012
        @test info.extras[:hopkins_per_dataset] == [0.81, 0.85]
    end

    @testset "cluster_statistics() abstract fallback errors" begin
        cam = IdealCamera(1:8, 1:8, 0.1)
        smld = BasicSMLD(SMLMData.Emitter2DFit{Float64}[], cam, 1, 1, Dict{String,Any}())
        # Fallback should refuse for any AbstractStatisticsConfig subtype without
        # a concrete cluster_statistics method, and the error message should
        # name HopkinsConfig as the available concrete backend.
        err = try
            cluster_statistics(smld, _DummyStatsCfg())
            nothing
        catch e
            e
        end
        @test err isa ErrorException
        @test occursin("HopkinsConfig", err.msg)
    end

end
