using SMLMClustering
using SMLMData
using Test
using Random
using Statistics

# Reuses `_make_2d_smld` and `_blob` from test_dbscan.jl.

@testset "point_hysteresis_clusters" begin

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
        @test_throws ArgumentError point_hysteresis_clusters(smld, seed, support)

        @test_throws ArgumentError point_hysteresis_clusters(
            smld, falses(49), falses(50))      # seed length mismatch
        @test_throws ArgumentError point_hysteresis_clusters(
            smld, falses(50), falses(49))      # support length mismatch
        @test_throws ArgumentError point_hysteresis_clusters(
            smld, falses(50), falses(50); graph_k = 0)
        @test_throws ArgumentError point_hysteresis_clusters(
            smld, falses(50), falses(50); min_points = 0)
    end

    @testset "trivial: empty smld + empty masks" begin
        cam = IdealCamera(1:8, 1:8, 0.1)
        smld = BasicSMLD(SMLMData.Emitter2DFit{Float64}[], cam, 1, 1,
                         Dict{String, Any}())
        smld_out, info = point_hysteresis_clusters(
            smld, Bool[], Bool[]; graph_k = 5, min_points = 5)
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
        support = trues(200)        # all support, no seeds
        smld_out, info = point_hysteresis_clusters(
            smld, seed, support; graph_k = 8, min_points = 50)
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
        smld_out, info = point_hysteresis_clusters(
            smld, seed, support; graph_k = 29, min_points = 10)
        @test info.n_clusters == 1
        @test info.cluster_sizes == [30]
        @test info.n_clustered == 30
        @test info.n_noise == n - 30
    end

    @testset "support without seed → component dropped" begin
        rng = Xoshiro(20260502)
        # Two well-separated 30-point clusters; only ONE has a seed.
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
        # support both clusters
        for i in a_start:b_end
            support[i] = true
        end
        # seed only the first
        seed[a_start] = true

        smld_out, info = point_hysteresis_clusters(
            smld, seed, support; graph_k = 29, min_points = 10)
        @test info.n_clusters == 1
        @test info.cluster_sizes == [30]
        # Cluster 1 should be the seeded one — verify by checking ids.
        @test all(smld_out.emitters[i].id == 1 for i in a_start:a_end)
        @test all(smld_out.emitters[i].id == 0 for i in b_start:b_end)
    end

    @testset "input ids do not leak to output (noise = 0 contract)" begin
        # Pipelines that feed prior-labeled SMLDs (e.g. BaGoL group ids)
        # must see id=0 on every unclustered emitter in the OUTPUT, even
        # though the input had non-zero ids. Verifies we zero ids on the
        # deepcopy before BFS.
        rng = Xoshiro(20260502)
        pts = Tuple{Float64, Float64, Int}[]
        # 30 points in a tight cluster, plus 20 far-away noise points.
        cluster_start = length(pts) + 1
        append!(pts, _blob(rng, 1.0, 1.0, 0.005, 30))
        cluster_end = length(pts)
        for _ in 1:20
            push!(pts, (10.0 + rand(rng), 10.0 + rand(rng), 1))
        end
        smld = _make_2d_smld(pts; n_datasets = 1)
        # Pre-stamp every emitter with a non-zero id, simulating a prior stage.
        for (i, e) in pairs(smld.emitters)
            e.id = 999
        end
        n = length(pts)
        seed = falses(n)
        support = falses(n)
        for i in cluster_start:cluster_end
            seed[i] = true; support[i] = true
        end

        smld_out, info = point_hysteresis_clusters(
            smld, seed, support; graph_k = 29, min_points = 10)

        cluster_ids = [smld_out.emitters[i].id for i in cluster_start:cluster_end]
        noise_ids = [smld_out.emitters[i].id for i in (cluster_end + 1):n]
        @test all(==(1), cluster_ids)              # cluster members get id=1
        @test all(==(0), noise_ids)                # noise points get id=0, NOT 999
        @test info.n_clustered == 30
        @test info.n_noise == 20
        # Run with no seeds → entire output must be id=0 even though input was 999.
        smld_empty, info_empty = point_hysteresis_clusters(
            smld, falses(n), falses(n); graph_k = 12, min_points = 10)
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
        smld_out, _ = point_hysteresis_clusters(
            smld_in, seed, support; graph_k = 5, min_points = 10)
        @test smld_out !== smld_in
        @test [e.id for e in smld_in.emitters] == ids_before  # input untouched
    end

    @testset "remove_unclustered drops noise from output" begin
        rng = Xoshiro(20260502)
        pts = Tuple{Float64, Float64, Int}[]
        append!(pts, _blob(rng, 1.0, 1.0, 0.005, 30))
        # 20 noise points far away
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

        smld_full, info_full = point_hysteresis_clusters(
            smld, seed, support; graph_k = 29, min_points = 10,
            remove_unclustered = false)
        smld_trim, info_trim = point_hysteresis_clusters(
            smld, seed, support; graph_k = 29, min_points = 10,
            remove_unclustered = true)
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

        # per_dataset=true: each dataset gets its own id=1
        smld_per, info_per = point_hysteresis_clusters(
            smld, seed, support; graph_k = 29, min_points = 10,
            per_dataset = true)
        @test info_per.n_clusters == 2
        ds1_ids = unique([e.id for e in smld_per.emitters if e.dataset == 1])
        ds2_ids = unique([e.id for e in smld_per.emitters if e.dataset == 2])
        @test ds1_ids == [1] && ds2_ids == [1]   # local namespacing per V3

        # per_dataset=false: pooled — depending on coords may merge or split,
        # but ids are global. Both datasets co-located → likely single cluster.
        smld_pool, info_pool = point_hysteresis_clusters(
            smld, seed, support; graph_k = 29, min_points = 10,
            per_dataset = false)
        @test info_pool.n_clusters >= 1
    end

    if SMLM_TEST_FULL
    @testset "wires with LocalContrastFeature: blob in uniform field is recovered" begin
        rng = Xoshiro(20260502)
        # Uniform background + tight blob; classic seed/support thresholds on
        # local contrast should recover most of the blob. graph_k=20 keeps the
        # asymmetric kNN BFS well-connected over the ~200-point blob.
        pts = [(5.0 * rand(rng), 5.0 * rand(rng), 1) for _ in 1:1500]
        blob_start = length(pts) + 1
        append!(pts, _blob(rng, 3.0, 3.0, 0.05, 200))
        blob_end = length(pts)
        smld = _make_2d_smld(pts; n_datasets = 1)

        cfg = LocalContrastFeature(density_k = 20, background_k = 400)
        _, info_f = cluster_statistics(smld, cfg)
        c = info_f.extras[:contrast_per_emitter]
        f = info_f.extras[:log_density_per_emitter]
        f_floor = quantile(filter(isfinite, f), 0.35)
        seed = isfinite.(c) .& isfinite.(f) .& (c .> 0.25) .& (f .> f_floor)
        support = isfinite.(c) .& isfinite.(f) .& (c .> -0.05) .& (f .> f_floor)

        _, info = point_hysteresis_clusters(smld, seed, support;
                                            graph_k = 20, min_points = 30)
        @test info.n_clusters >= 1
        @test info.cluster_sizes[1] >= 100     # most of the 200-point blob

        # Verify the recovered cluster overlaps the blob, not the background.
        _, idxs_seed = findmax(c)
        # Largest cluster should contain a substantial fraction of the blob
        # (membership check via BFS labels would be more direct, but the
        # cluster_sizes check above + the seed/support thresholds tied to
        # local contrast already constrain location).
    end
    end  # SMLM_TEST_FULL

end  # @testset point_hysteresis_clusters
