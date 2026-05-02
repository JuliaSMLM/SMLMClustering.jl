using SMLMClustering
using SMLMData
using Test

# Test-tier gate. Default off → only the fast tier runs (interface, exports,
# Config round-trip, basic correctness). Set SMLM_TEST_FULL=true (or any value
# accepted as truthy here) to also run the thorough tier — multi-blob ground-
# truth recovery, edge cases, statistical value checks, large-n stress, etc.
# Cross-package convention shared with SMLMAnalysis / SMLMBaGoL /
# SMLMDriftCorrection; see Round History in STATUS.md for the rollout note.
const SMLM_TEST_FULL = lowercase(get(ENV, "SMLM_TEST_FULL", "false")) in ("true", "1", "yes")

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
    include("test_hdbscan.jl")
    include("test_hierarchical.jl")
    include("test_voronoi.jl")
    include("test_cluster_statistics.jl")
    include("test_hopkins.jl")
    include("test_voronoi_density.jl")
    include("test_local_contrast.jl")
    include("test_mrf_density.jl")
    include("test_point_hysteresis.jl")

end

if !SMLM_TEST_FULL
    @info "Skipping thorough tests; set SMLM_TEST_FULL=true to enable"
end
