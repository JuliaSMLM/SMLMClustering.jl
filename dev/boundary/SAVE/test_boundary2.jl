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
            2,                                    # n_clusters
            [(1, 1), (1, 2)],                     # cluster_keys
            [1.0, 2.5],                           # area
            [4.0, 6.3],                           # perimeter
            [[1, 2, 3, 1], [1, 2, 3, 1]],          # indices
            [[1, 2, 3, 1], [1, 2, 3, 1]],          # indicesConvex
            [[1, 2, 3, 1], [1, 2, 3, 1]],          # indicesTight
            [1.0, 2.5],                           # areaConvex
            [0.8, 2.0],                           # areaTight
            [4.0, 6.3],                           # perimeterConvex
            [3.6, 5.8],                           # perimeterTight
            [(0.5, 0.5), (2.0, 2.0)],              # center
            [sqrt(1.0 / pi) / 2, sqrt(2.5 / pi) / 2],  # equiv_radius
            [0.7, 0.8],                           # compactness
            [0.6, 0.75],                          # circularity
            [4.0 / 3.6, 6.3 / 5.8],                # convexity
            [0.8, 0.8],                           # solidity
            [5, 8],                               # n_pts_per_cluster
            [5.0, 3.2],                           # n_pts_per_area
            [0.3, 0.5],                           # sigma_actual
            [1.2, 2.1],                           # cluster_width
            [3.0, 3.0],                           # min_c2c_dists
            [2.0, 2.0],                           # min_e2e_dists
            3.0,                                  # min_c2c_dist
            2.0,                                  # min_e2e_dist
            [[0.1, 0.1, 0.1, 0.1, 0.1], fill(0.2, 8)],  # nn_within_clusters
            0.001,                                 # elapsed_s
        )
        @test sinfo.n_clusters == 2
        @test sinfo.cluster_keys == [(1, 1), (1, 2)]
        @test sinfo.area == [1.0, 2.5]
        @test sinfo.perimeter == [4.0, 6.3]
        @test sinfo.indicesConvex == [[1, 2, 3, 1], [1, 2, 3, 1]]
        @test sinfo.center == [(0.5, 0.5), (2.0, 2.0)]
        @test sinfo.n_pts_per_cluster == [5, 8]
        @test sinfo.min_c2c_dist ≈ 3.0
        @test sinfo.min_e2e_dist ≈ 2.0
        @test length(sinfo.nn_within_clusters[2]) == 8
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
        @test isempty(sinfo.indices)
        @test isempty(sinfo.center)
        @test isnan(sinfo.min_c2c_dist)
        @test isnan(sinfo.min_e2e_dist)
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
        # smld has no emitters, so this cluster key has zero matching raw
        # points — every raw-point-derived field degrades to NaN/empty
        # rather than throwing.
        @test sinfo.n_pts_per_cluster[1] == 0
        @test sinfo.indices[1] == Int[]
        @test isnan(sinfo.sigma_actual[1])
        @test isempty(sinfo.nn_within_clusters[1])
        @test sinfo.center[1] === (NaN, NaN) || all(isnan, sinfo.center[1])

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
        # equiv_radius/compactness derive from area/perimeter alone, so they
        # are well-defined even though smld carries no matching raw points.
        @test sinfo.equiv_radius[1] ≈ sqrt(1.0 / pi) / 2 atol = 1e-12
        @test sinfo.compactness[1] ≈ pi / 4 atol = 1e-12
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

        # Real matching points → raw-point-derived stats are well-defined.
        @test sinfo.n_pts_per_cluster == [6]
        @test all(x -> x isa Real && isfinite(x), sinfo.center[1])
        @test all(0.0 .< sinfo.compactness .<= 1.0)
        @test all(0.0 .< sinfo.circularity .<= 1.0)
        @test all(0.0 .< sinfo.solidity .<= 1.0)
        # convexity = perimeterConvex / perimeterTight is not bounded above by
        # 1 in general: the shrink=1 "tightest" alpha shape keeps only the
        # largest boundary loop, which can trace a small subset of the
        # cluster's points (here a 3-point sub-triangle) rather than an
        # enclosing curve around every point, making its perimeter shorter
        # than the convex hull's.
        @test all(x -> isfinite(x) && x > 0, sinfo.convexity)
        @test length(sinfo.nn_within_clusters[1]) == 6
        # A single cluster has no "other cluster" to measure against.
        @test isnan(sinfo.min_c2c_dist)
        @test isnan(sinfo.min_e2e_dist)
    end
    end  # SMLM_TEST_FULL

end
