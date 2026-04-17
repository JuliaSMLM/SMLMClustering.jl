using SMLMClustering
using SMLMData
using Test
using Random

# Reuse `_make_2d_smld` and `_blob` from test_dbscan.jl (included earlier
# in runtests.jl so those helpers are already defined at top-level).

@testset "Voronoi backend" begin

    @testset "config construction" begin
        cfg = VoronoiConfig()
        @test cfg isa AbstractClusterConfig
        @test cfg.density_factor == 2.0
        @test cfg.min_points == 5
        @test cfg.use_3d === false
        @test cfg.per_dataset === true
        @test cfg.remove_unclustered === false

        cfg2 = VoronoiConfig(density_factor = 3.5, min_points = 3,
                              per_dataset = false, remove_unclustered = true)
        @test cfg2.density_factor == 3.5
        @test cfg2.min_points == 3
        @test cfg2.per_dataset === false
        @test cfg2.remove_unclustered === true
    end

    @testset "three well-separated blobs + scattered noise" begin
        rng = Xoshiro(20260417)
        σ = 0.010             # 10 nm — tight clusters
        n_per_blob = 60       # more points per blob → cleaner density statistic

        pts = Tuple{Float64,Float64,Int}[]
        append!(pts, _blob(rng, 2.0, 2.0, σ, n_per_blob))
        append!(pts, _blob(rng, 4.0, 2.0, σ, n_per_blob))
        append!(pts, _blob(rng, 3.0, 4.0, σ, n_per_blob))
        # 60 scattered "noise" points spread over a 6×6 μm region — their
        # cells are much larger than blob cells, so they are not "dense".
        for _ in 1:60
            push!(pts, (6.0 * rand(rng), 6.0 * rand(rng), 1))
        end
        smld = _make_2d_smld(pts; n_datasets = 1)
        n_in = length(smld.emitters)

        cfg = VoronoiConfig(density_factor = 2.0, min_points = 5,
                             per_dataset = false)
        smld_out, info = cluster(smld, cfg)

        @test info isa ClusterInfo
        @test info.algorithm === :voronoi
        @test info.n_locs_in == n_in
        @test info.elapsed_s >= 0
        @test info.n_clusters == 3
        @test sum(info.cluster_sizes) == info.n_clustered
        @test info.n_clustered + info.n_noise == info.n_locs_in
        # Blob points have cell areas ≪ mean; the great majority should
        # cluster. Some blob-edge cells clipped by nearby noise points can
        # exceed the density threshold, so allow a small loss margin.
        @test info.n_clustered >= 3 * n_per_blob - 20
        # Noise points should be largely excluded.
        @test info.n_noise >= 40
    end

    @testset "labels written to emitter.id + remove_unclustered" begin
        rng = Xoshiro(11)
        pts = Tuple{Float64,Float64,Int}[]
        append!(pts, _blob(rng, 1.0, 1.0, 0.005, 40))
        append!(pts, _blob(rng, 5.0, 5.0, 0.005, 40))
        # 20 scattered points, area cells much larger → noise.
        for _ in 1:20
            push!(pts, (10.0 * rand(rng), 10.0 * rand(rng), 1))
        end
        smld = _make_2d_smld(pts; n_datasets = 1)
        n_in = length(smld.emitters)

        cfg = VoronoiConfig(density_factor = 2.0, min_points = 5,
                             per_dataset = false)
        smld_keep, info_keep = cluster(smld, cfg)

        @test length(smld_keep.emitters) == n_in
        @test info_keep.n_clusters == 2
        @test all(e -> e.id in 0:info_keep.n_clusters, smld.emitters)
        @test any(e -> e.id == 0, smld.emitters)
        for k in 1:info_keep.n_clusters
            @test count(e -> e.id == k, smld.emitters) == info_keep.cluster_sizes[k]
        end

        smld2 = _make_2d_smld(pts; n_datasets = 1)
        cfg_rm = VoronoiConfig(density_factor = 2.0, min_points = 5,
                                per_dataset = false, remove_unclustered = true)
        smld_rm, info_rm = cluster(smld2, cfg_rm)
        @test info_rm.n_clusters == 2
        @test length(smld_rm.emitters) == info_rm.n_clustered
        @test all(e -> e.id != 0, smld_rm.emitters)
    end

    @testset "per_dataset label namespace is local" begin
        rng = Xoshiro(21)
        pts = Tuple{Float64,Float64,Int}[]
        append!(pts, _blob(rng, 1.0, 1.0, 0.005, 30; dataset = 1))
        append!(pts, _blob(rng, 3.0, 3.0, 0.005, 30; dataset = 1))
        # scattered per-dataset noise so blobs stand out in the density map
        for _ in 1:20
            push!(pts, (5.0 * rand(rng), 5.0 * rand(rng), 1))
        end
        append!(pts, _blob(rng, 1.0, 1.0, 0.005, 30; dataset = 2))
        append!(pts, _blob(rng, 3.0, 3.0, 0.005, 30; dataset = 2))
        for _ in 1:20
            push!(pts, (5.0 * rand(rng), 5.0 * rand(rng), 2))
        end
        smld = _make_2d_smld(pts; n_datasets = 2)

        cfg = VoronoiConfig(density_factor = 2.0, min_points = 5,
                             per_dataset = true)
        _, info = cluster(smld, cfg)
        @test info.n_clusters == 4

        ids_ds1 = sort!(unique(e.id for e in smld.emitters if e.dataset == 1 && e.id > 0))
        ids_ds2 = sort!(unique(e.id for e in smld.emitters if e.dataset == 2 && e.id > 0))
        @test ids_ds1 == [1, 2]
        @test ids_ds2 == [1, 2]

        # Contrast: per_dataset=false merges same-coordinate blobs across datasets.
        smld_flat = _make_2d_smld(pts; n_datasets = 2)
        cfg_flat = VoronoiConfig(density_factor = 2.0, min_points = 5,
                                   per_dataset = false)
        _, info_flat = cluster(smld_flat, cfg_flat)
        @test info_flat.n_clusters == 2
    end

    @testset "argument validation" begin
        smld = _make_2d_smld([(0.0, 0.0, 1), (1.0, 0.0, 1), (0.0, 1.0, 1)])
        @test_throws ArgumentError cluster(smld, VoronoiConfig(density_factor = 0.0))
        @test_throws ArgumentError cluster(smld, VoronoiConfig(density_factor = -1.0))
        @test_throws ArgumentError cluster(smld, VoronoiConfig(min_points = 0))
        @test_throws ArgumentError cluster(smld, VoronoiConfig(use_3d = true))
    end

    @testset "duplicate coordinates raise ArgumentError" begin
        # Exact-coincident (x,y) generators cause DelaunayTriangulation.get_area
        # to raise KeyError; we want a clean ArgumentError before triangulation.
        pts = [(0.0, 0.0, 1), (1.0, 0.0, 1), (0.0, 1.0, 1), (1.0, 1.0, 1),
               (0.0, 0.0, 1)]  # duplicate of first point
        smld = _make_2d_smld(pts; n_datasets = 1)
        cfg = VoronoiConfig(density_factor = 2.0, min_points = 1, per_dataset = false)
        @test_throws ArgumentError cluster(smld, cfg)

        # Duplicate in a multi-dataset split: per_dataset=true means each dataset
        # is processed independently, so the error fires per-group.
        pts2 = [(0.0, 0.0, 1), (1.0, 0.0, 1), (0.0, 1.0, 1), (0.0, 0.0, 1)]
        smld2 = _make_2d_smld(pts2; n_datasets = 1)
        cfg2 = VoronoiConfig(density_factor = 2.0, min_points = 1, per_dataset = true)
        @test_throws ArgumentError cluster(smld2, cfg2)
    end

    @testset "degenerate groups (<3 points) are all noise" begin
        # One dataset with only 2 points; Voronoi can't tessellate → all noise.
        smld = _make_2d_smld([(0.0, 0.0, 1), (1.0, 1.0, 1)]; n_datasets = 1)
        cfg = VoronoiConfig(density_factor = 2.0, min_points = 1, per_dataset = false)
        smld_out, info = cluster(smld, cfg)
        @test info.n_clusters == 0
        @test info.n_clustered == 0
        @test info.n_noise == 2
        @test all(e -> e.id == 0, smld_out.emitters)
    end

    @testset "empty SMLD" begin
        smld = _make_2d_smld(Tuple{Float64,Float64,Int}[]; n_datasets = 1)
        cfg = VoronoiConfig(per_dataset = false)
        smld_out, info = cluster(smld, cfg)
        @test info.n_locs_in == 0
        @test info.n_clusters == 0
        @test info.n_clustered == 0
        @test info.n_noise == 0
        @test isempty(info.cluster_sizes)
        @test isempty(smld_out.emitters)
    end

    @testset "density_factor threshold behavior" begin
        # Tight blob (area ≪ mean) + sparse halo. With density_factor=2, blob
        # survives. With density_factor = a very small fraction, the threshold
        # becomes huge and almost every point passes as "dense" → the blob
        # absorbs the halo into a single large cluster.
        rng = Xoshiro(42)
        pts = Tuple{Float64,Float64,Int}[]
        append!(pts, _blob(rng, 2.0, 2.0, 0.005, 50))
        for _ in 1:50
            push!(pts, (4.0 * rand(rng), 4.0 * rand(rng), 1))
        end
        n_in = length(pts)

        smld_strict = _make_2d_smld(pts; n_datasets = 1)
        cfg_strict = VoronoiConfig(density_factor = 2.0, min_points = 5,
                                    per_dataset = false)
        _, info_strict = cluster(smld_strict, cfg_strict)
        @test info_strict.n_clusters == 1
        @test info_strict.n_clustered < n_in   # halo points remain noise

        smld_loose = _make_2d_smld(pts; n_datasets = 1)
        cfg_loose = VoronoiConfig(density_factor = 0.01, min_points = 5,
                                   per_dataset = false)
        _, info_loose = cluster(smld_loose, cfg_loose)
        # Extremely permissive threshold → (almost) every point is dense, so
        # one giant component forms; more points get clustered than under strict.
        @test info_loose.n_clustered > info_strict.n_clustered
    end

end
