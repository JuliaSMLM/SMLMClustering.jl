using SMLMClustering
using SMLMData
using Test
using Random

# Helper: build a 2D BasicSMLD of Emitter2DFit{Float64} from (x, y, dataset)
# tuples (microns). frame=1 throughout.
function _make_2d_smld_hdbscan(points::Vector{Tuple{Float64,Float64,Int}};
                              n_datasets::Int = maximum(p[3] for p in points))
    cam = IdealCamera(1:64, 1:64, 0.1)
    emitters = [Emitter2DFit{Float64}(
        x, y, 1000.0, 10.0, 0.01, 0.01, 50.0, 2.0;
        frame = 1, dataset = ds,
    ) for (x, y, ds) in points]
    BasicSMLD(emitters, cam, 1, n_datasets, Dict{String,Any}())
end

function _blob_hdbscan(rng, cx, cy, σ, n; dataset = 1)
    [(cx + σ * randn(rng), cy + σ * randn(rng), dataset) for _ in 1:n]
end

@testset "HDBSCAN backend" begin

    @testset "config construction" begin
        cfg = HDBSCANConfig()
        @test cfg isa AbstractClusterConfig
        @test cfg.min_points == 5
        @test cfg.min_cluster_size === nothing
        @test cfg.knn_graph_k == 30
        @test cfg.cluster_selection_method === :eom
        @test cfg.allow_single_cluster === false
        @test cfg.use_3d === false
        @test cfg.per_dataset === true
        @test cfg.remove_unclustered === false

        cfg2 = HDBSCANConfig(min_points = 10, min_cluster_size = 20,
                             knn_graph_k = 50, cluster_selection_method = :leaf,
                             allow_single_cluster = true,
                             use_3d = true, per_dataset = false,
                             remove_unclustered = true)
        @test cfg2.min_points == 10
        @test cfg2.min_cluster_size == 20
        @test cfg2.knn_graph_k == 50
        @test cfg2.cluster_selection_method === :leaf
        @test cfg2.allow_single_cluster === true
        @test cfg2.use_3d === true
        @test cfg2.per_dataset === false
        @test cfg2.remove_unclustered === true
    end

    @testset "argument validation" begin
        smld = _make_2d_smld_hdbscan([(0.0, 0.0, 1)])
        @test_throws ArgumentError cluster(smld, HDBSCANConfig(min_points = 0))
        @test_throws ArgumentError cluster(smld, HDBSCANConfig(knn_graph_k = 0))
        @test_throws ArgumentError cluster(smld, HDBSCANConfig(min_cluster_size = 1))
        @test_throws ArgumentError cluster(smld, HDBSCANConfig(cluster_selection_method = :bogus))
    end

    if SMLM_TEST_FULL
    @testset "two well-separated Gaussian blobs of different density" begin
        rng = Xoshiro(20260428)
        # Blob A: 60 points, σ=15 nm at (1, 1) μm.
        # Blob B: 80 points, σ= 8 nm at (3, 1) μm  (denser).
        # Centers 2 μm apart — vastly larger than either σ. Real split at root.
        pts = Tuple{Float64,Float64,Int}[]
        append!(pts, _blob_hdbscan(rng, 1.0, 1.0, 0.015, 60))
        append!(pts, _blob_hdbscan(rng, 3.0, 1.0, 0.008, 80))
        smld = _make_2d_smld_hdbscan(pts; n_datasets = 1)

        cfg = HDBSCANConfig(min_points = 5, min_cluster_size = 15,
                            knn_graph_k = 30, per_dataset = false)
        smld_out, info = cluster(smld, cfg)
        @test info isa ClusterInfo
        @test info.algorithm === :hdbscan
        @test info.n_clusters == 2
        # Both clusters should pull most of their points
        @test info.cluster_sizes[1] >= 50
        @test info.cluster_sizes[2] >= 50
        # Persistence reported via metadata
        pers = smld_out.metadata["hdbscan_cluster_persistence"]
        @test length(pers) == 2
        @test all(p -> p >= 0, pers)
        @test all(isfinite, pers)
    end
    end  # SMLM_TEST_FULL

    if SMLM_TEST_FULL
    @testset "two blobs + tight satellite below mcs (hex-vs-structure)" begin
        # Two well-separated structure-scale blobs (real split at root); a tight
        # 6-point satellite attached to one of them must be treated as noise
        # (below min_cluster_size). EOM should keep the two parent blobs and
        # drop the satellite.
        rng = Xoshiro(20260428 + 1)
        pts = Tuple{Float64,Float64,Int}[]
        append!(pts, _blob_hdbscan(rng, 0.0, 0.0, 0.020, 100))
        append!(pts, _blob_hdbscan(rng, 2.0, 0.0, 0.020, 100))
        # 6-point hex-scale satellite far from both — should be noise (mcs=15).
        append!(pts, _blob_hdbscan(rng, 4.0, 0.0, 0.001, 6))
        smld = _make_2d_smld_hdbscan(pts; n_datasets = 1)

        cfg = HDBSCANConfig(min_points = 5, min_cluster_size = 15,
                            knn_graph_k = 30, per_dataset = false)
        smld_out, info = cluster(smld, cfg)
        @test info.n_clusters == 2
        # Each big blob keeps most of its 100 points.
        @test all(s -> s >= 80, info.cluster_sizes)
        # The 6 satellite points are noise.
        @test info.n_noise >= 6
        pers = smld_out.metadata["hdbscan_cluster_persistence"]
        @test length(pers) == 2
        @test all(p -> p >= 0, pers)
    end
    end  # SMLM_TEST_FULL

    if SMLM_TEST_FULL
    @testset "uniform random — sparse means few/no clusters" begin
        rng = Xoshiro(20260428 + 2)
        # 200 points uniform in [0, 5]² μm. Sparse, no real clusters. With
        # default allow_single_cluster=false the EOM rule cannot pick root,
        # so we expect ZERO selected clusters (everything noise).
        pts = Tuple{Float64,Float64,Int}[]
        for _ in 1:200
            push!(pts, (5.0 * rand(rng), 5.0 * rand(rng), 1))
        end
        smld = _make_2d_smld_hdbscan(pts; n_datasets = 1)

        cfg = HDBSCANConfig(min_points = 5, min_cluster_size = 30,
                            knn_graph_k = 30, per_dataset = false)
        _, info = cluster(smld, cfg)
        @test info.n_clusters <= 3
        @test info.n_clustered + info.n_noise == info.n_locs_in
    end
    end  # SMLM_TEST_FULL

    if SMLM_TEST_FULL
    @testset "persistence is non-negative for every selected cluster" begin
        rng = Xoshiro(20260428 + 3)
        pts = Tuple{Float64,Float64,Int}[]
        # Three planted clusters with different σ — separated enough to give
        # real splits.
        append!(pts, _blob_hdbscan(rng, 1.0, 1.0, 0.010, 50))
        append!(pts, _blob_hdbscan(rng, 3.0, 1.0, 0.020, 50))
        append!(pts, _blob_hdbscan(rng, 1.0, 3.0, 0.030, 50))
        smld = _make_2d_smld_hdbscan(pts; n_datasets = 1)

        cfg = HDBSCANConfig(min_points = 5, min_cluster_size = 20,
                            knn_graph_k = 30, per_dataset = false)
        smld_out, info = cluster(smld, cfg)
        @test info.n_clusters >= 1
        pers = smld_out.metadata["hdbscan_cluster_persistence"]
        births = smld_out.metadata["hdbscan_cluster_lambda_birth"]
        @test length(pers) == info.n_clusters
        @test length(births) == info.n_clusters
        @test all(p -> p >= 0, pers)
        @test all(isfinite, pers)
        @test all(b -> b >= 0, births)
    end
    end  # SMLM_TEST_FULL

    if SMLM_TEST_FULL
    @testset "labels written to emitter.id (non-mutating input)" begin
        rng = Xoshiro(20260428 + 4)
        pts = Tuple{Float64,Float64,Int}[]
        append!(pts, _blob_hdbscan(rng, 1.0, 1.0, 0.005, 40))
        append!(pts, _blob_hdbscan(rng, 5.0, 5.0, 0.005, 40))
        smld = _make_2d_smld_hdbscan(pts; n_datasets = 1)

        cfg = HDBSCANConfig(min_points = 5, min_cluster_size = 15,
                            knn_graph_k = 30, per_dataset = false)
        smld_out, info = cluster(smld, cfg)
        @test info.n_clusters == 2
        @test all(e -> e.id == 0, smld.emitters)  # input untouched
        @test all(e -> e.id in 0:info.n_clusters, smld_out.emitters)
        for k in 1:info.n_clusters
            @test count(e -> e.id == k, smld_out.emitters) == info.cluster_sizes[k]
        end
    end
    end  # SMLM_TEST_FULL

    if SMLM_TEST_FULL
    @testset "remove_unclustered drops noise" begin
        # Two blobs with extra noise outliers (below mcs as a group) → 2 clusters,
        # outliers dropped on output.
        rng = Xoshiro(20260428 + 5)
        pts = Tuple{Float64,Float64,Int}[]
        append!(pts, _blob_hdbscan(rng, 0.0, 0.0, 0.015, 50))
        append!(pts, _blob_hdbscan(rng, 3.0, 0.0, 0.015, 50))
        # 5 distant outliers
        for k in 1:5
            push!(pts, (10.0 + k, 10.0 + k, 1))
        end
        smld = _make_2d_smld_hdbscan(pts; n_datasets = 1)

        cfg = HDBSCANConfig(min_points = 5, min_cluster_size = 15,
                            knn_graph_k = 30, per_dataset = false,
                            remove_unclustered = true)
        smld_out, info = cluster(smld, cfg)
        @test info.n_clusters == 2
        @test length(smld_out.emitters) == info.n_clustered
        @test all(e -> e.id != 0, smld_out.emitters)
    end
    end  # SMLM_TEST_FULL

    if SMLM_TEST_FULL
    @testset "allow_single_cluster=true returns one cluster on a tight blob" begin
        # Single tight blob (no internal real splits). Default EOM returns
        # zero clusters; allow_single_cluster=true returns 1.
        rng = Xoshiro(20260428 + 6)
        pts = _blob_hdbscan(rng, 1.0, 1.0, 0.005, 60)
        smld = _make_2d_smld_hdbscan(pts; n_datasets = 1)

        cfg_default = HDBSCANConfig(min_points = 5, min_cluster_size = 15,
                                    knn_graph_k = 30, per_dataset = false)
        _, info_def = cluster(smld, cfg_default)
        @test info_def.n_clusters == 0

        cfg_single = HDBSCANConfig(min_points = 5, min_cluster_size = 15,
                                   knn_graph_k = 30, per_dataset = false,
                                   allow_single_cluster = true)
        _, info_sin = cluster(smld, cfg_single)
        @test info_sin.n_clusters == 1
        @test info_sin.cluster_sizes[1] == 60
    end
    end  # SMLM_TEST_FULL

    if SMLM_TEST_FULL
    @testset "kNN graph small still completes via component-bridging" begin
        # Two well-separated blobs with knn_graph_k=2. The kNN graph alone
        # cannot bridge the blobs, but the bridging step in the backend
        # finds the cheapest cross-component edge and continues. The end
        # result should still be 2 clusters.
        rng = Xoshiro(20260428 + 7)
        pts = Tuple{Float64,Float64,Int}[]
        append!(pts, _blob_hdbscan(rng, 0.0, 0.0, 0.005, 60))
        append!(pts, _blob_hdbscan(rng, 5.0, 5.0, 0.005, 60))
        smld = _make_2d_smld_hdbscan(pts; n_datasets = 1)

        cfg = HDBSCANConfig(min_points = 5, min_cluster_size = 15,
                            knn_graph_k = 2, per_dataset = false)
        _, info = cluster(smld, cfg)
        @test info.n_clusters == 2
    end
    end  # SMLM_TEST_FULL

end
