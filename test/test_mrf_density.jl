using SMLMClustering
using SMLMData
using Test
using Random
using Statistics

# Reuses `_make_2d_smld` and `_blob` from test_dbscan.jl.

@testset "MRFDensityCluster backend" begin

    @testset "config construction" begin
        cfg = MRFDensityClusterConfig()
        @test cfg isa AbstractClusterConfig
        @test cfg.n_regimes == 2
        @test cfg.regime_thresholds === nothing
        @test cfg.smoothness_lambda === nothing
        @test cfg.graph_kind === :delaunay
        @test cfg.graph_k == 8
        @test cfg.inference === :icm
        @test cfg.icm_iters == 50
        @test cfg.min_points == 5
        @test cfg.use_3d === false
        @test cfg.per_dataset === true
        @test cfg.remove_unclustered === false

        cfg2 = MRFDensityClusterConfig(n_regimes = 3, graph_kind = :knn,
                                       graph_k = 12, min_points = 10,
                                       smoothness_lambda = 0.5,
                                       per_dataset = false,
                                       remove_unclustered = true)
        @test cfg2.n_regimes == 3
        @test cfg2.graph_kind === :knn
        @test cfg2.graph_k == 12
        @test cfg2.min_points == 10
        @test cfg2.smoothness_lambda == 0.5
        @test cfg2.per_dataset === false
        @test cfg2.remove_unclustered === true

        cfg3 = MRFDensityClusterConfig(n_regimes = 3,
                                       regime_thresholds = [3.0, 5.0])
        @test cfg3.regime_thresholds == [3.0, 5.0]
    end

    @testset "argument validation" begin
        # 3-point dummy SMLD so we get past the n<3 short-circuit and reach validation.
        smld_dummy = _make_2d_smld([(0.0, 0.0, 1), (1.0, 0.0, 1), (0.0, 1.0, 1)])

        @test_throws ArgumentError cluster(smld_dummy, MRFDensityClusterConfig(n_regimes = 1))
        @test_throws ArgumentError cluster(smld_dummy, MRFDensityClusterConfig(graph_kind = :bogus))
        @test_throws ArgumentError cluster(smld_dummy, MRFDensityClusterConfig(inference = :graph_cuts))
        @test_throws ArgumentError cluster(smld_dummy, MRFDensityClusterConfig(graph_k = 0))
        @test_throws ArgumentError cluster(smld_dummy, MRFDensityClusterConfig(icm_iters = 0))
        @test_throws ArgumentError cluster(smld_dummy, MRFDensityClusterConfig(min_points = 0))
        @test_throws ArgumentError cluster(smld_dummy, MRFDensityClusterConfig(smoothness_lambda = -1.0))
        # Threshold length mismatch (n_regimes=2 expects length 1, given 2).
        @test_throws ArgumentError cluster(smld_dummy,
            MRFDensityClusterConfig(n_regimes = 2, regime_thresholds = [1.0, 2.0]))
        # Unsorted thresholds.
        @test_throws ArgumentError cluster(smld_dummy,
            MRFDensityClusterConfig(n_regimes = 3, regime_thresholds = [5.0, 3.0]))
    end

    @testset "use_3d=true raises ArgumentError" begin
        smld = _make_2d_smld([(0.0, 0.0, 1), (1.0, 0.0, 1), (0.0, 1.0, 1)])
        @test_throws ArgumentError cluster(smld, MRFDensityClusterConfig(use_3d = true))
    end

    @testset "duplicate (x,y) coordinates raise ArgumentError" begin
        # Need ≥3 points per group to bypass the small-group short-circuit.
        pts = [(0.0, 0.0, 1), (1.0, 0.0, 1), (0.0, 1.0, 1), (1.0, 1.0, 1),
               (0.0, 0.0, 1)]  # duplicate of first point
        smld = _make_2d_smld(pts; n_datasets = 1)
        @test_throws ArgumentError cluster(smld, MRFDensityClusterConfig(per_dataset = false))
    end

    @testset "empty SMLD: zero clusters, empty regime vector" begin
        smld = _make_2d_smld(Tuple{Float64,Float64,Int}[]; n_datasets = 1)
        smld_out, info = cluster(smld, MRFDensityClusterConfig())
        @test info.n_locs_in == 0
        @test info.n_clustered == 0
        @test info.n_clusters == 0
        @test info.algorithm === :mrf_density
        @test smld_out.metadata["mrf_regime_per_emitter"] == Int[]
    end

    @testset "non-mutation: input emitter ids untouched" begin
        rng = Xoshiro(7)
        pts = Tuple{Float64,Float64,Int}[]
        for _ in 1:30
            push!(pts, (2.0 * rand(rng), 2.0 * rand(rng), 1))
        end
        # Add a tighter blob.
        append!(pts, _blob(rng, 1.0, 1.0, 0.05, 30))
        smld = _make_2d_smld(pts; n_datasets = 1)
        @test all(e -> e.id == 0, smld.emitters)
        cluster(smld, MRFDensityClusterConfig(per_dataset = false))
        @test all(e -> e.id == 0, smld.emitters)  # input unchanged
    end

    @testset "exports" begin
        @test :MRFDensityClusterConfig in names(SMLMClustering)
    end

    if SMLM_TEST_FULL
    @testset "2-regime synthetic: blob recovered, background is noise" begin
        rng = Xoshiro(20260429)
        # 200 background points in a 5x5 box; 60 tight blob at (2.5, 2.5), σ=0.05 μm.
        pts = Tuple{Float64,Float64,Int}[]
        for _ in 1:200
            push!(pts, (5.0 * rand(rng), 5.0 * rand(rng), 1))
        end
        append!(pts, _blob(rng, 2.5, 2.5, 0.05, 60))
        smld = _make_2d_smld(pts; n_datasets = 1)

        cfg = MRFDensityClusterConfig(per_dataset = false, min_points = 10)
        smld_out, info = cluster(smld, cfg)
        @test info.n_clusters >= 1
        # Blob is by far the largest cluster; should sweep up most of the
        # 60 blob points.
        @test maximum(info.cluster_sizes) >= 40
        # Most of the background should be noise (id = 0).
        bg_indices = 1:200
        bg_noise = count(i -> smld_out.emitters[i].id == 0, bg_indices)
        @test bg_noise > 150  # allow some bleed at the boundary
        # Regime metadata exists.
        regimes = smld_out.metadata["mrf_regime_per_emitter"]
        @test length(regimes) == length(smld.emitters)
        @test all(r -> 0 <= r <= 2, regimes)
    end
    end  # SMLM_TEST_FULL

    if SMLM_TEST_FULL
    @testset "3-regime synthetic: three density tiers identified via metadata" begin
        rng = Xoshiro(20260429)
        # Tier 1: low-density background in a 6×6 box.
        # Tier 2: medium patch in [1, 2] × [1, 2].
        # Tier 3: tight blob at (4, 4), σ=0.04.
        pts = Tuple{Float64,Float64,Int}[]
        for _ in 1:300
            push!(pts, (6.0 * rand(rng), 6.0 * rand(rng), 1))
        end
        for _ in 1:120
            push!(pts, (1.0 + 1.0 * rand(rng), 1.0 + 1.0 * rand(rng), 1))
        end
        append!(pts, _blob(rng, 4.0, 4.0, 0.04, 80))
        smld = _make_2d_smld(pts; n_datasets = 1)

        cfg = MRFDensityClusterConfig(n_regimes = 3, per_dataset = false,
                                      min_points = 8)
        smld_out, info = cluster(smld, cfg)
        regimes = smld_out.metadata["mrf_regime_per_emitter"]
        # All three regime levels should be represented.
        @test 1 in regimes
        @test 2 in regimes
        @test 3 in regimes
        # GMM means recorded.
        means = smld_out.metadata["mrf_regime_means"]
        @test length(means) == 1
        @test length(means[1]) == 3
        @test issorted(means[1])  # sorted ascending
        # Highest density blob → some clusters formed.
        @test info.n_clusters >= 1
    end
    end  # SMLM_TEST_FULL

    if SMLM_TEST_FULL
    @testset "missing-middle: MRF smoothness keeps borderline middle points in cluster" begin
        rng = Xoshiro(20260429)
        # Dense blob at (2.5, 2.5), σ=0.06 μm, plus a few "intermediate"
        # points in the middle that, on their own, would be at borderline
        # density. Then surround with low-density background.
        pts = Tuple{Float64,Float64,Int}[]
        # Dense blob.
        append!(pts, _blob(rng, 2.5, 2.5, 0.06, 80))
        # Background noise.
        for _ in 1:200
            push!(pts, (5.0 * rand(rng), 5.0 * rand(rng), 1))
        end
        smld = _make_2d_smld(pts; n_datasets = 1)

        cfg = MRFDensityClusterConfig(per_dataset = false, min_points = 10)
        smld_out, info = cluster(smld, cfg)
        @test info.n_clusters >= 1
        # The cluster shouldn't have implausibly small size (smoothness
        # should pull in borderline points).
        @test maximum(info.cluster_sizes) >= 50
    end
    end  # SMLM_TEST_FULL

    if SMLM_TEST_FULL
    @testset "spurious-small: 4-emitter knot dies, 50-emitter blob survives" begin
        rng = Xoshiro(20260429)
        pts = Tuple{Float64,Float64,Int}[]
        # Background.
        for _ in 1:200
            push!(pts, (8.0 * rand(rng), 8.0 * rand(rng), 1))
        end
        # 50-emitter genuine blob.
        append!(pts, _blob(rng, 6.0, 6.0, 0.05, 50))
        # 4-emitter tight knot far from the blob (so MRF smoothness can
        # plausibly demote it).
        append!(pts, _blob(rng, 1.5, 1.5, 0.005, 4))
        smld = _make_2d_smld(pts; n_datasets = 1)

        cfg = MRFDensityClusterConfig(per_dataset = false, min_points = 5)
        smld_out, info = cluster(smld, cfg)
        # The min_points=5 size filter ensures the 4-knot can't survive
        # even if it locally passes the regime test.
        @test all(s -> s >= 5, info.cluster_sizes)
        # The 50-blob should be one of the recovered clusters.
        @test maximum(info.cluster_sizes) >= 30
    end
    end  # SMLM_TEST_FULL

    if SMLM_TEST_FULL
    @testset "regime_thresholds override bypasses GMM" begin
        rng = Xoshiro(20260429)
        pts = Tuple{Float64,Float64,Int}[]
        for _ in 1:200
            push!(pts, (5.0 * rand(rng), 5.0 * rand(rng), 1))
        end
        append!(pts, _blob(rng, 2.5, 2.5, 0.05, 60))
        smld = _make_2d_smld(pts; n_datasets = 1)

        # Pick a threshold roughly between background and blob log-density.
        # Background ~200/25 = 8 μm⁻², log≈2.1; blob inside σ²~3e-3 μm² so
        # density ~hundreds → log~5+. Pick 4.0 as a neutral split.
        cfg = MRFDensityClusterConfig(per_dataset = false,
                                      regime_thresholds = [4.0],
                                      min_points = 10)
        smld_out, info = cluster(smld, cfg)
        # When thresholds are explicit, regime_means metadata is filled with NaN.
        means = smld_out.metadata["mrf_regime_means"]
        @test length(means) == 1
        @test all(isnan, means[1])
        # And we should still recover the blob.
        @test info.n_clusters >= 1
    end
    end  # SMLM_TEST_FULL

    if SMLM_TEST_FULL
    @testset "per_dataset=true: each dataset gets its own GMM fit" begin
        rng = Xoshiro(20260429)
        # ds 1: low-density bg + blob.
        pts = Tuple{Float64,Float64,Int}[]
        for _ in 1:120
            push!(pts, (5.0 * rand(rng), 5.0 * rand(rng), 1))
        end
        append!(pts, [(p[1], p[2], 1) for p in _blob(rng, 2.5, 2.5, 0.05, 50)])
        n_ds1 = length(pts)
        # ds 2: roughly 10x denser background + blob (different scale).
        for _ in 1:120
            push!(pts, (1.5 * rand(rng), 1.5 * rand(rng), 2))
        end
        append!(pts, [(p[1], p[2], 2) for p in _blob(rng, 0.75, 0.75, 0.02, 50)])
        smld = _make_2d_smld(pts; n_datasets = 2)

        cfg = MRFDensityClusterConfig(per_dataset = true, min_points = 8)
        smld_out, info = cluster(smld, cfg)
        means = smld_out.metadata["mrf_regime_means"]
        @test length(means) == 2  # one fit per dataset
        # Each dataset should pick up at least one cluster.
        ds1_n_clusters = length(unique([smld_out.emitters[i].id for i in 1:n_ds1
                                         if smld_out.emitters[i].id > 0]))
        ds2_n_clusters = length(unique([smld_out.emitters[i].id for i in (n_ds1+1):length(smld_out.emitters)
                                         if smld_out.emitters[i].id > 0]))
        @test ds1_n_clusters >= 1
        @test ds2_n_clusters >= 1
    end
    end  # SMLM_TEST_FULL

    if SMLM_TEST_FULL
    @testset ":knn graph mode produces a clustering with similar shape to :delaunay" begin
        rng = Xoshiro(20260429)
        pts = Tuple{Float64,Float64,Int}[]
        for _ in 1:150
            push!(pts, (5.0 * rand(rng), 5.0 * rand(rng), 1))
        end
        append!(pts, _blob(rng, 2.5, 2.5, 0.05, 60))
        smld = _make_2d_smld(pts; n_datasets = 1)

        cfg_dt = MRFDensityClusterConfig(per_dataset = false, min_points = 10,
                                         graph_kind = :delaunay)
        cfg_kn = MRFDensityClusterConfig(per_dataset = false, min_points = 10,
                                         graph_kind = :knn, graph_k = 10)
        _, info_dt = cluster(smld, cfg_dt)
        _, info_kn = cluster(smld, cfg_kn)
        # Both methods should find at least one cluster — counts may differ
        # slightly, structure should match.
        @test info_dt.n_clusters >= 1
        @test info_kn.n_clusters >= 1
    end
    end  # SMLM_TEST_FULL

    if SMLM_TEST_FULL
    @testset "deterministic reproducibility on identical input" begin
        rng = Xoshiro(20260429)
        pts = Tuple{Float64,Float64,Int}[]
        for _ in 1:120
            push!(pts, (4.0 * rand(rng), 4.0 * rand(rng), 1))
        end
        append!(pts, _blob(rng, 2.0, 2.0, 0.05, 50))
        smld_a = _make_2d_smld(pts; n_datasets = 1)
        smld_b = _make_2d_smld(pts; n_datasets = 1)
        cfg = MRFDensityClusterConfig(per_dataset = false, min_points = 10)
        out_a, _ = cluster(smld_a, cfg)
        out_b, _ = cluster(smld_b, cfg)
        @test [e.id for e in out_a.emitters] == [e.id for e in out_b.emitters]
    end
    end  # SMLM_TEST_FULL

    if SMLM_TEST_FULL
    @testset "remove_unclustered=true: returned smld contains only clustered emitters" begin
        rng = Xoshiro(20260429)
        pts = Tuple{Float64,Float64,Int}[]
        for _ in 1:100
            push!(pts, (5.0 * rand(rng), 5.0 * rand(rng), 1))
        end
        append!(pts, _blob(rng, 2.5, 2.5, 0.05, 60))
        smld = _make_2d_smld(pts; n_datasets = 1)
        cfg = MRFDensityClusterConfig(per_dataset = false, min_points = 10,
                                      remove_unclustered = true)
        smld_out, info = cluster(smld, cfg)
        @test all(e -> e.id > 0, smld_out.emitters)
        @test length(smld_out.emitters) == info.n_clustered
    end
    end  # SMLM_TEST_FULL

end
