using SMLMClustering
using SMLMData
using Test
using Random

# Reuse helpers from test_dbscan.jl (loaded first in runtests.jl).
# _make_2d_smld and _blob are defined there.

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

    if SMLM_TEST_FULL
    @testset "three well-separated blobs + noise (single linkage)" begin
        rng = Xoshiro(20260417)
        σ = 0.010             # 10 nm — tight clusters
        n_per_blob = 40

        pts = Tuple{Float64,Float64,Int}[]
        append!(pts, _blob(rng, 2.0, 2.0, σ, n_per_blob))
        append!(pts, _blob(rng, 4.0, 2.0, σ, n_per_blob))
        append!(pts, _blob(rng, 3.0, 4.0, σ, n_per_blob))
        for _ in 1:30
            push!(pts, (6.0 * rand(rng), 6.0 * rand(rng), 1))
        end
        smld = _make_2d_smld(pts; n_datasets = 1)
        n_in = length(smld.emitters)

        # Single linkage at 100 nm: blobs stay separate (nearest-neighbor within
        # blob ≪ 100 nm; nearest inter-blob distance ≫ 100 nm).
        cfg = HierarchicalConfig(cut_threshold = 100.0, linkage = :single,
                                  min_points = 5, per_dataset = false)
        smld_out, info = cluster(smld, cfg)

        @test info isa ClusterInfo
        @test info.algorithm === :hierarchical
        @test info.n_locs_in == n_in
        @test info.elapsed_s >= 0
        @test info.n_clusters == 3
        @test sum(info.cluster_sizes) == info.n_clustered
        @test info.n_clustered + info.n_noise == info.n_locs_in
        # All 120 blob points should cluster; the 30 scattered noise points
        # form singletons (< min_points=5) and are tagged noise.
        @test info.n_clustered >= 115
        @test info.n_noise >= 20
    end
    end  # SMLM_TEST_FULL

    if SMLM_TEST_FULL
    @testset "n_clusters path (Ward linkage)" begin
        rng = Xoshiro(20260418)
        σ = 0.010
        n_per_blob = 40

        pts = Tuple{Float64,Float64,Int}[]
        append!(pts, _blob(rng, 2.0, 2.0, σ, n_per_blob))
        append!(pts, _blob(rng, 4.0, 2.0, σ, n_per_blob))
        append!(pts, _blob(rng, 3.0, 4.0, σ, n_per_blob))
        smld = _make_2d_smld(pts; n_datasets = 1)

        # Ward into 3 clusters: each blob is a cluster.
        cfg = HierarchicalConfig(n_clusters = 3, linkage = :ward,
                                  min_points = 5, per_dataset = false)
        smld_out, info = cluster(smld, cfg)
        @test info.n_clusters == 3
        @test info.n_clustered == 3 * n_per_blob
        @test info.n_noise == 0
    end
    end  # SMLM_TEST_FULL

    if SMLM_TEST_FULL
    @testset "labels written to emitter.id + remove_unclustered + non-mutating" begin
        rng = Xoshiro(10)
        pts = Tuple{Float64,Float64,Int}[]
        append!(pts, _blob(rng, 1.0, 1.0, 0.005, 30))
        append!(pts, _blob(rng, 5.0, 5.0, 0.005, 30))
        for k in 1:5
            push!(pts, (10.0 + k, 10.0 + k, 1))
        end
        smld = _make_2d_smld(pts; n_datasets = 1)
        n_in = length(smld.emitters)

        cfg = HierarchicalConfig(cut_threshold = 100.0, linkage = :single,
                                  min_points = 5, per_dataset = false)
        smld_keep, info_keep = cluster(smld, cfg)

        @test length(smld_keep.emitters) == n_in
        @test info_keep.n_clusters == 2
        # Input SMLD emitters are NOT mutated (non-mutating semantics).
        @test all(e -> e.id == 0, smld.emitters)
        # Output SMLD carries the labels.
        @test all(e -> e.id in 0:info_keep.n_clusters, smld_keep.emitters)
        @test any(e -> e.id == 0, smld_keep.emitters)
        for k in 1:info_keep.n_clusters
            @test count(e -> e.id == k, smld_keep.emitters) == info_keep.cluster_sizes[k]
        end

        cfg_rm = HierarchicalConfig(cut_threshold = 100.0, linkage = :single,
                                     min_points = 5, per_dataset = false,
                                     remove_unclustered = true)
        smld_rm, info_rm = cluster(smld, cfg_rm)
        @test info_rm.n_clusters == 2
        @test length(smld_rm.emitters) == info_rm.n_clustered
        @test all(e -> e.id != 0, smld_rm.emitters)
        # Original smld still untouched after second call.
        @test all(e -> e.id == 0, smld.emitters)
    end
    end  # SMLM_TEST_FULL

    if SMLM_TEST_FULL
    @testset "per_dataset label namespace is local" begin
        rng = Xoshiro(20)
        pts = Tuple{Float64,Float64,Int}[]
        append!(pts, _blob(rng, 1.0, 1.0, 0.005, 20; dataset = 1))
        append!(pts, _blob(rng, 3.0, 3.0, 0.005, 20; dataset = 1))
        append!(pts, _blob(rng, 1.0, 1.0, 0.005, 20; dataset = 2))
        append!(pts, _blob(rng, 3.0, 3.0, 0.005, 20; dataset = 2))
        smld = _make_2d_smld(pts; n_datasets = 2)

        cfg = HierarchicalConfig(cut_threshold = 100.0, linkage = :single,
                                  min_points = 5, per_dataset = true)
        smld_out, info = cluster(smld, cfg)
        @test info.n_clusters == 4
        @test info.n_clustered == 80

        ids_ds1 = sort!(unique(e.id for e in smld_out.emitters if e.dataset == 1))
        ids_ds2 = sort!(unique(e.id for e in smld_out.emitters if e.dataset == 2))
        @test ids_ds1 == [1, 2]
        @test ids_ds2 == [1, 2]

        cfg_flat = HierarchicalConfig(cut_threshold = 100.0, linkage = :single,
                                       min_points = 5, per_dataset = false)
        _, info_flat = cluster(smld, cfg_flat)
        @test info_flat.n_clusters == 2
        @test info_flat.n_clustered == 80
    end
    end  # SMLM_TEST_FULL

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

    if SMLM_TEST_FULL
    @testset "3D clustering path" begin
        cam = IdealCamera(1:64, 1:64, 0.1)
        rng = Xoshiro(30)
        emitters = SMLMData.Emitter3DFit{Float64}[]
        for _ in 1:25
            push!(emitters, Emitter3DFit{Float64}(
                2.0 + 0.005 * randn(rng), 2.0 + 0.005 * randn(rng), 0.0 + 0.005 * randn(rng),
                1000.0, 10.0, 0.01, 0.01, 0.02, 50.0, 2.0; frame = 1, dataset = 1))
        end
        for _ in 1:25
            push!(emitters, Emitter3DFit{Float64}(
                2.0 + 0.005 * randn(rng), 2.0 + 0.005 * randn(rng), 1.0 + 0.005 * randn(rng),
                1000.0, 10.0, 0.01, 0.01, 0.02, 50.0, 2.0; frame = 1, dataset = 1))
        end
        smld3 = BasicSMLD(emitters, cam, 1, 1, Dict{String,Any}())

        cfg = HierarchicalConfig(cut_threshold = 100.0, linkage = :single,
                                  min_points = 5, use_3d = true, per_dataset = false)
        _, info = cluster(smld3, cfg)
        @test info.n_clusters == 2
        @test info.n_clustered == 50
    end
    end  # SMLM_TEST_FULL

    if SMLM_TEST_FULL
    @testset "empty SMLD" begin
        smld = _make_2d_smld(Tuple{Float64,Float64,Int}[]; n_datasets = 1)
        cfg = HierarchicalConfig(cut_threshold = 100.0, per_dataset = false)
        smld_out, info = cluster(smld, cfg)
        @test info.n_locs_in == 0
        @test info.n_clusters == 0
        @test info.n_clustered == 0
        @test info.n_noise == 0
        @test isempty(info.cluster_sizes)
        @test isempty(smld_out.emitters)
    end
    end  # SMLM_TEST_FULL

    if SMLM_TEST_FULL
    @testset "min_points filters small clusters as noise" begin
        # Three tight blobs of 10 pts each + two singletons.
        rng = Xoshiro(40)
        pts = Tuple{Float64,Float64,Int}[]
        for (cx, cy) in [(1.0, 1.0), (3.0, 1.0), (2.0, 3.0)]
            append!(pts, _blob(rng, cx, cy, 0.005, 10))
        end
        push!(pts, (10.0, 10.0, 1), (11.0, 11.0, 1))
        smld = _make_2d_smld(pts; n_datasets = 1)

        # With min_points=5: blobs (10 pts each) survive, singletons → noise.
        cfg5 = HierarchicalConfig(cut_threshold = 50.0, linkage = :single,
                                   min_points = 5, per_dataset = false)
        _, info5 = cluster(smld, cfg5)
        @test info5.n_clusters == 3
        @test info5.n_noise == 2

        # With min_points=15: all blobs (10 pts each) are below threshold → all noise.
        cfg15 = HierarchicalConfig(cut_threshold = 50.0, linkage = :single,
                                    min_points = 15, per_dataset = false)
        _, info15 = cluster(smld, cfg15)
        @test info15.n_clusters == 0
        @test info15.n_noise == length(pts)
    end
    end  # SMLM_TEST_FULL

end
