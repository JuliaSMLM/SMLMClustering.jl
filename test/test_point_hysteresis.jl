using SMLMClustering
using SMLMData
using Test
using Random
using Statistics

# Reuses `_make_2d_smld` and `_blob` from test_dbscan.jl.

@testset "PointHysteresisConfig backend" begin

    @testset "config construction" begin
        cfg = PointHysteresisConfig()
        @test cfg isa AbstractClusterConfig
        @test cfg.graph_k == 12
        @test cfg.min_points == 100
        @test cfg.use_3d === false
        @test cfg.per_dataset === false
        @test cfg.remove_unclustered === false

        cfg2 = PointHysteresisConfig(graph_k = 20, min_points = 50,
                                     per_dataset = true)
        @test cfg2.graph_k == 20
        @test cfg2.min_points == 50
        @test cfg2.per_dataset === true
    end

    @testset "argument validation" begin
        rng = Xoshiro(20260502)
        pts = [(0.1 * randn(rng), 0.1 * randn(rng), 1) for _ in 1:50]
        smld = _make_2d_smld(pts; n_datasets = 1)
        seed = falses(50)
        support = falses(50)
        seed[1] = true
        support[1] = true
        seed[2] = true
        # support[2] left false → seed without support: should error
        @test_throws ArgumentError cluster(smld, PointHysteresisConfig();
                                            seed = seed, support = support)

        @test_throws ArgumentError cluster(smld, PointHysteresisConfig();
                                            seed = falses(49), support = falses(50))
        @test_throws ArgumentError cluster(smld, PointHysteresisConfig();
                                            seed = falses(50), support = falses(49))
        @test_throws ArgumentError cluster(smld,
            PointHysteresisConfig(graph_k = 0);
            seed = falses(50), support = falses(50))
        @test_throws ArgumentError cluster(smld,
            PointHysteresisConfig(min_points = 0);
            seed = falses(50), support = falses(50))
    end

    @testset "trivial: empty smld + empty masks" begin
        cam = IdealCamera(1:8, 1:8, 0.1)
        smld = BasicSMLD(SMLMData.Emitter2DFit{Float64}[], cam, 1, 1,
                         Dict{String, Any}())
        smld_out, info = cluster(smld, PointHysteresisConfig(graph_k = 5,
                                                              min_points = 5);
                                  seed = Bool[], support = Bool[])
        @test info.n_locs_in == 0
        @test info.n_clusters == 0
        @test info.algorithm === :point_hysteresis
        @test length(smld_out.emitters) == 0
    end

    @testset "no seeds → no clusters" begin
        rng = Xoshiro(20260502)
        pts = [(0.5 * randn(rng), 0.5 * randn(rng), 1) for _ in 1:200]
        smld = _make_2d_smld(pts; n_datasets = 1)
        seed = falses(200)
        support = trues(200)
        smld_out, info = cluster(smld,
            PointHysteresisConfig(graph_k = 8, min_points = 50);
            seed = seed, support = support)
        @test info.n_clusters == 0
        @test info.n_clustered == 0
        @test info.n_noise == 200
        @test all(e -> e.id == 0, smld_out.emitters)
    end

    @testset "min_points filter rejects small components" begin
        rng = Xoshiro(20260502)
        # Two well-separated clusters, sizes 30 and 5, plus noise.
        # graph_k=29 gives a full-mesh graph on the 30-point cluster so BFS
        # is guaranteed to reach every support point in it (the asymmetric
        # k-NN BFS can leave a couple of outlier points unreachable for
        # tight blobs when k < n-1; this is correct behavior of the
        # algorithm and matches the dev-script reference).
        pts = Tuple{Float64, Float64, Int}[]
        big_start = length(pts) + 1
        append!(pts, _blob(rng, 1.0, 1.0, 0.005, 30))
        big_end = length(pts)
        small_start = length(pts) + 1
        append!(pts, _blob(rng, 4.0, 4.0, 0.005, 5))
        small_end = length(pts)

        smld = _make_2d_smld(pts; n_datasets = 1)
        n = length(pts)
        seed = falses(n)
        support = falses(n)
        for i in big_start:big_end
            seed[i] = true; support[i] = true
        end
        for i in small_start:small_end
            seed[i] = true; support[i] = true
        end
        smld_out, info = cluster(smld,
            PointHysteresisConfig(graph_k = 29, min_points = 10);
            seed = seed, support = support)
        @test info.n_clusters == 1
        @test info.cluster_sizes == [30]
        @test info.n_clustered == 30
        @test info.n_noise == n - 30
    end

    @testset "support without seed → component dropped" begin
        rng = Xoshiro(20260502)
        pts = Tuple{Float64, Float64, Int}[]
        a_start = length(pts) + 1
        append!(pts, _blob(rng, 1.0, 1.0, 0.005, 30))
        a_end = length(pts)
        b_start = length(pts) + 1
        append!(pts, _blob(rng, 4.0, 4.0, 0.005, 30))
        b_end = length(pts)

        smld = _make_2d_smld(pts; n_datasets = 1)
        n = length(pts)
        seed = falses(n)
        support = falses(n)
        for i in a_start:b_end
            support[i] = true
        end
        seed[a_start] = true

        smld_out, info = cluster(smld,
            PointHysteresisConfig(graph_k = 29, min_points = 10);
            seed = seed, support = support)
        @test info.n_clusters == 1
        @test info.cluster_sizes == [30]
        @test all(smld_out.emitters[i].id == 1 for i in a_start:a_end)
        @test all(smld_out.emitters[i].id == 0 for i in b_start:b_end)
    end

    @testset "input ids do not leak to output (noise = 0 contract)" begin
        rng = Xoshiro(20260502)
        pts = Tuple{Float64, Float64, Int}[]
        cluster_start = length(pts) + 1
        append!(pts, _blob(rng, 1.0, 1.0, 0.005, 30))
        cluster_end = length(pts)
        for _ in 1:20
            push!(pts, (10.0 + rand(rng), 10.0 + rand(rng), 1))
        end
        smld = _make_2d_smld(pts; n_datasets = 1)
        for (i, e) in pairs(smld.emitters)
            e.id = 999
        end
        n = length(pts)
        seed = falses(n)
        support = falses(n)
        for i in cluster_start:cluster_end
            seed[i] = true; support[i] = true
        end

        smld_out, info = cluster(smld,
            PointHysteresisConfig(graph_k = 29, min_points = 10);
            seed = seed, support = support)

        cluster_ids = [smld_out.emitters[i].id for i in cluster_start:cluster_end]
        noise_ids = [smld_out.emitters[i].id for i in (cluster_end + 1):n]
        @test all(==(1), cluster_ids)
        @test all(==(0), noise_ids)
        @test info.n_clustered == 30
        @test info.n_noise == 20

        smld_empty, info_empty = cluster(smld,
            PointHysteresisConfig(graph_k = 12, min_points = 10);
            seed = falses(n), support = falses(n))
        @test all(e -> e.id == 0, smld_empty.emitters)
        @test info_empty.n_clustered == 0
    end

    @testset "non-mutation of input" begin
        rng = Xoshiro(20260502)
        pts = Tuple{Float64, Float64, Int}[]
        append!(pts, _blob(rng, 1.0, 1.0, 0.005, 30))
        smld_in = _make_2d_smld(pts; n_datasets = 1)
        ids_before = [e.id for e in smld_in.emitters]
        seed = trues(length(pts))
        support = trues(length(pts))
        smld_out, _ = cluster(smld_in,
            PointHysteresisConfig(graph_k = 5, min_points = 10);
            seed = seed, support = support)
        @test smld_out !== smld_in
        @test [e.id for e in smld_in.emitters] == ids_before
    end

    @testset "remove_unclustered drops noise from output" begin
        rng = Xoshiro(20260502)
        pts = Tuple{Float64, Float64, Int}[]
        append!(pts, _blob(rng, 1.0, 1.0, 0.005, 30))
        for _ in 1:20
            push!(pts, (10.0 + rand(rng), 10.0 + rand(rng), 1))
        end
        smld = _make_2d_smld(pts; n_datasets = 1)
        n = length(pts)
        seed = falses(n)
        support = falses(n)
        for i in 1:30
            seed[i] = true; support[i] = true
        end

        smld_full, info_full = cluster(smld,
            PointHysteresisConfig(graph_k = 29, min_points = 10,
                                  remove_unclustered = false);
            seed = seed, support = support)
        smld_trim, info_trim = cluster(smld,
            PointHysteresisConfig(graph_k = 29, min_points = 10,
                                  remove_unclustered = true);
            seed = seed, support = support)
        @test length(smld_full.emitters) == n
        @test length(smld_trim.emitters) == 30
        @test info_full.n_clusters == 1 == info_trim.n_clusters
    end

    @testset "per_dataset gives independent label namespaces" begin
        rng = Xoshiro(20260502)
        pts = Tuple{Float64, Float64, Int}[]
        append!(pts, _blob(rng, 1.0, 1.0, 0.005, 30; dataset = 1))
        append!(pts, _blob(rng, 1.0, 1.0, 0.005, 30; dataset = 2))
        smld = _make_2d_smld(pts; n_datasets = 2)
        n = length(pts)
        seed = trues(n)
        support = trues(n)

        smld_per, info_per = cluster(smld,
            PointHysteresisConfig(graph_k = 29, min_points = 10,
                                  per_dataset = true);
            seed = seed, support = support)
        @test info_per.n_clusters == 2
        ds1_ids = unique([e.id for e in smld_per.emitters if e.dataset == 1])
        ds2_ids = unique([e.id for e in smld_per.emitters if e.dataset == 2])
        @test ds1_ids == [1] && ds2_ids == [1]

        smld_pool, info_pool = cluster(smld,
            PointHysteresisConfig(graph_k = 29, min_points = 10,
                                  per_dataset = false);
            seed = seed, support = support)
        @test info_pool.n_clusters >= 1
    end

    @testset "config is reusable across multiple SMLDs" begin
        # Verify the "config = stable reusable knobs" invariant: same cfg
        # instance, two different SMLDs with their own seed/support, both
        # produce a result.
        rng = Xoshiro(20260502)
        cfg = PointHysteresisConfig(graph_k = 29, min_points = 10)

        pts1 = _blob(rng, 1.0, 1.0, 0.005, 30)
        smld1 = _make_2d_smld(pts1; n_datasets = 1)
        s1 = trues(30); sup1 = trues(30)

        pts2 = _blob(rng, 5.0, 5.0, 0.005, 30)
        smld2 = _make_2d_smld(pts2; n_datasets = 1)
        s2 = trues(30); sup2 = trues(30)

        _, info1 = cluster(smld1, cfg; seed = s1, support = sup1)
        _, info2 = cluster(smld2, cfg; seed = s2, support = sup2)
        @test info1.n_clusters == 1
        @test info2.n_clusters == 1
    end

    if SMLM_TEST_FULL
    @testset "wires with LocalContrastFeature: blob in uniform field is recovered" begin
        rng = Xoshiro(20260502)
        pts = [(5.0 * rand(rng), 5.0 * rand(rng), 1) for _ in 1:1500]
        blob_start = length(pts) + 1
        append!(pts, _blob(rng, 3.0, 3.0, 0.05, 200))
        blob_end = length(pts)
        smld = _make_2d_smld(pts; n_datasets = 1)

        feat_cfg = LocalContrastFeature(density_k = 20, background_k = 400)
        _, info_f = cluster_statistics(smld, feat_cfg)
        c = info_f.extras[:contrast_per_emitter]
        f = info_f.extras[:log_density_per_emitter]
        f_floor = quantile(filter(isfinite, f), 0.35)
        seed = isfinite.(c) .& isfinite.(f) .& (c .> 0.25) .& (f .> f_floor)
        support = isfinite.(c) .& isfinite.(f) .& (c .> -0.05) .& (f .> f_floor)

        _, info = cluster(smld,
            PointHysteresisConfig(graph_k = 20, min_points = 30);
            seed = seed, support = support)
        @test info.n_clusters >= 1
        @test info.cluster_sizes[1] >= 100
    end
    end  # SMLM_TEST_FULL

end  # @testset PointHysteresisConfig backend
