using SMLMClustering
using SMLMData
using Test

# Six-point test cluster used throughout this file (coordinates in μm).
# Indices: 1:(0,0)  2:(1,1)  3:(2,1)  4:(3,0)  5:(3,2)  6:(0,2)
# Convex hull vertices are {1,4,5,6} = {(0,0),(3,0),(3,2),(0,2)}.
# Interior points are {2,3} = {(1,1),(2,1)}.
const _BND_PTS_X = [0.0, 1.0, 2.0, 3.0, 3.0, 0.0]
const _BND_PTS_Y = [0.0, 1.0, 1.0, 0.0, 2.0, 2.0]

@testset "boundary backend" begin

    @testset "config construction" begin
        cfg = BoundaryClustersConfig()
        @test cfg.shrink == 0.5
        @test cfg.per_dataset === true

        cfg2 = BoundaryClustersConfig(shrink = 0.0, per_dataset = false)
        @test cfg2.shrink == 0.0
        @test cfg2.per_dataset === false

        cfg3 = BoundaryClustersConfig(shrink = 1.0)
        @test cfg3.shrink == 1.0
    end

    @testset "ClusterBoundaryInfo construction" begin
        info = ClusterBoundaryInfo(
            2,
            [(1, 1), (1, 2)],
            [[0.0, 1.0, 0.0], [2.0, 3.0, 2.0]],
            [[0.0, 0.0, 0.0], [1.0, 2.0, 1.0]],
            0.5,
            0.001,
        )
        @test info.n_clusters == 2
        @test info.cluster_keys == [(1, 1), (1, 2)]
        @test info.shrink == 0.5
        @test info.elapsed_s ≈ 0.001
    end

    @testset "boundary_cluster edge cases" begin
        @test boundary_cluster(Float64[], Float64[]) == Int[]
        @test boundary_cluster([1.0], [1.0]) == [1]
        @test boundary_cluster([1.0, 2.0], [0.0, 0.0]) == [1, 2, 1]
        result3 = boundary_cluster([0.0, 1.0, 0.5], [0.0, 0.0, 1.0])
        @test result3 == [1, 2, 3, 1]
    end

    @testset "boundary_cluster argument errors" begin
        x3 = [0.0, 1.0, 0.5]
        y3 = [0.0, 0.0, 1.0]
        @test_throws ArgumentError boundary_cluster(x3, y3; shrink = -0.1)
        @test_throws ArgumentError boundary_cluster(x3, y3; shrink = 1.1)
        @test_throws ArgumentError boundary_cluster([1.0, 2.0], [1.0])
    end

    @testset "boundary_clusters argument errors" begin
        cam = IdealCamera(1:8, 1:8, 0.1)
        smld_empty = BasicSMLD(SMLMData.Emitter2DFit{Float64}[], cam, 1, 1,
                               Dict{String,Any}())
        @test_throws ArgumentError boundary_clusters(
            smld_empty, BoundaryClustersConfig(shrink = -0.1))
        @test_throws ArgumentError boundary_clusters(
            smld_empty, BoundaryClustersConfig(shrink = 1.5))
    end

    if SMLM_TEST_FULL
    @testset "boundary_cluster shrink=0 — convex hull" begin
        x = _BND_PTS_X
        y = _BND_PTS_Y
        result = boundary_cluster(x, y; shrink = 0.0)

        @test result isa Vector{Int}
        @test result[1] == result[end]                  # closed polygon
        unique_idx = unique(result)
        @test length(unique_idx) == 4                   # hull has 4 vertices
        # Interior points (1,1) and (2,1) must not appear on the convex hull.
        @test 2 ∉ unique_idx
        @test 3 ∉ unique_idx
        @test SMLMClustering._signed_area_polygon(x, y, result) > 0   # CCW
    end
    end  # SMLM_TEST_FULL

    if SMLM_TEST_FULL
    @testset "boundary_cluster shrink=0.5 — intermediate alpha shape" begin
        x = _BND_PTS_X
        y = _BND_PTS_Y
        result = boundary_cluster(x, y; shrink = 0.5)

        @test result isa Vector{Int}
        @test result[1] == result[end]
        @test length(unique(result)) in 4:6
        @test SMLMClustering._signed_area_polygon(x, y, result) > 0
    end
    end  # SMLM_TEST_FULL

    if SMLM_TEST_FULL
    @testset "boundary_cluster shrink=1 — tightest alpha shape" begin
        x = _BND_PTS_X
        y = _BND_PTS_Y
        result = boundary_cluster(x, y; shrink = 1.0)

        @test result isa Vector{Int}
        @test result[1] == result[end]
        @test length(result) >= 4
        @test SMLMClustering._signed_area_polygon(x, y, result) > 0
    end
    end  # SMLM_TEST_FULL

    if SMLM_TEST_FULL
    @testset "boundary_clusters — single cluster, three shrink values" begin
        pts = Tuple{Float64,Float64,Int}[
            (0.0, 0.0, 1), (1.0, 1.0, 1), (2.0, 1.0, 1),
            (3.0, 0.0, 1), (3.0, 2.0, 1), (0.0, 2.0, 1),
        ]
        smld = _make_2d_smld(pts; n_datasets = 1)
        # eps=5000 nm >> 3.6 μm diagonal → all 6 points form one cluster.
        smld_out, _ = cluster(smld,
            DBSCANConfig(eps_nm = 5000.0, min_points = 2, per_dataset = false))

        for shrink in (0.0, 0.5, 1.0)
            @testset "shrink=$shrink" begin
                cfg_bnd = BoundaryClustersConfig(shrink = shrink, per_dataset = true)
                (smld_ref, binfo) = boundary_clusters(smld_out, cfg_bnd)

                @test smld_ref === smld_out             # pass-through: same reference
                @test binfo isa ClusterBoundaryInfo
                @test binfo.n_clusters == 1
                @test binfo.cluster_keys == [(1, 1)]
                @test binfo.shrink == shrink
                @test binfo.elapsed_s >= 0

                bx = binfo.boundaries_x[1]
                by = binfo.boundaries_y[1]
                @test bx[1] == bx[end]                 # closed x
                @test by[1] == by[end]                 # closed y
                @test length(bx) >= 4

                if shrink == 0.0
                    # Convex hull: only extreme coordinates appear.
                    @test all(v -> v ∈ (0.0, 3.0), bx)
                    @test all(v -> v ∈ (0.0, 2.0), by)
                end
            end
        end
    end
    end  # SMLM_TEST_FULL

    @testset "ClusterBoundaryStatisticsInfo construction" begin
        sinfo = ClusterBoundaryStatisticsInfo(
            2,
            [(1, 1), (1, 2)],
            [1.0, 2.5],
            [4.0, 6.3],
            0.001,
        )
        @test sinfo.n_clusters == 2
        @test sinfo.cluster_keys == [(1, 1), (1, 2)]
        @test sinfo.area == [1.0, 2.5]
        @test sinfo.perimeter == [4.0, 6.3]
        @test sinfo.elapsed_s ≈ 0.001
    end

    @testset "boundary_statistics — pass-through SMLD and empty binfo" begin
        cam        = IdealCamera(1:8, 1:8, 0.1)
        smld_empty = BasicSMLD(SMLMData.Emitter2DFit{Float64}[], cam, 1, 1,
                               Dict{String,Any}())
        binfo_empty = ClusterBoundaryInfo(0, Tuple{Int,Int}[],
                                          Vector{Float64}[], Vector{Float64}[],
                                          0.5, 0.0)
        (smld_ref, sinfo) = boundary_statistics(smld_empty, binfo_empty)
        @test smld_ref === smld_empty           # pass-through: same reference
        @test sinfo isa ClusterBoundaryStatisticsInfo
        @test sinfo.n_clusters == 0
        @test isempty(sinfo.area)
        @test isempty(sinfo.perimeter)
        @test sinfo.elapsed_s >= 0
    end

    @testset "boundary_statistics — degenerate short polygon → NaN" begin
        cam  = IdealCamera(1:8, 1:8, 0.1)
        smld = BasicSMLD(SMLMData.Emitter2DFit{Float64}[], cam, 1, 1,
                         Dict{String,Any}())
        # Single-point polygon: length == 1 (< 4)
        binfo = ClusterBoundaryInfo(1, [(1, 1)],
                                    [[0.5]], [[0.5]],
                                    0.5, 0.0)
        (_, sinfo) = boundary_statistics(smld, binfo)
        @test isnan(sinfo.area[1])
        @test isnan(sinfo.perimeter[1])

        # Two-point (degenerate segment): [1,2,1] has length 3 (< 4)
        binfo2 = ClusterBoundaryInfo(1, [(1, 1)],
                                     [[0.0, 1.0, 0.0]], [[0.0, 0.0, 0.0]],
                                     0.5, 0.0)
        (_, sinfo2) = boundary_statistics(smld, binfo2)
        @test isnan(sinfo2.area[1])
        @test isnan(sinfo2.perimeter[1])
    end

    if SMLM_TEST_FULL
    @testset "boundary_statistics — unit square polygon" begin
        cam  = IdealCamera(1:8, 1:8, 0.1)
        smld = BasicSMLD(SMLMData.Emitter2DFit{Float64}[], cam, 1, 1,
                         Dict{String,Any}())
        # Unit square (CCW): (0,0)→(1,0)→(1,1)→(0,1)→(0,0)
        bx = [0.0, 1.0, 1.0, 0.0, 0.0]
        by = [0.0, 0.0, 1.0, 1.0, 0.0]
        binfo = ClusterBoundaryInfo(1, [(1, 1)], [bx], [by], 0.5, 0.0)
        (_, sinfo) = boundary_statistics(smld, binfo)

        @test sinfo.n_clusters == 1
        @test sinfo.cluster_keys == [(1, 1)]
        @test sinfo.area[1] ≈ 1.0 atol = 1e-12
        @test sinfo.perimeter[1] ≈ 4.0 atol = 1e-12
    end
    end  # SMLM_TEST_FULL

    if SMLM_TEST_FULL
    @testset "boundary_statistics — integration: cluster → boundary_clusters → boundary_statistics" begin
        pts = Tuple{Float64,Float64,Int}[
            (0.0, 0.0, 1), (1.0, 1.0, 1), (2.0, 1.0, 1),
            (3.0, 0.0, 1), (3.0, 2.0, 1), (0.0, 2.0, 1),
        ]
        smld = _make_2d_smld(pts; n_datasets = 1)
        smld_out, _ = cluster(smld,
            DBSCANConfig(eps_nm = 5000.0, min_points = 2, per_dataset = false))
        _, binfo = boundary_clusters(smld_out, BoundaryClustersConfig(shrink = 0.0))
        smld_ref, sinfo = boundary_statistics(smld_out, binfo)

        @test smld_ref === smld_out             # pass-through
        @test sinfo isa ClusterBoundaryStatisticsInfo
        @test sinfo.n_clusters == binfo.n_clusters
        @test sinfo.cluster_keys == binfo.cluster_keys
        @test all(a -> a > 0, sinfo.area)
        @test all(p -> p > 0, sinfo.perimeter)
        @test sinfo.elapsed_s >= 0
    end
    end  # SMLM_TEST_FULL

end
