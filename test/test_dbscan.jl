using SMLMClustering
using SMLMData
using Test
using Random

# Helper: build a BasicSMLD of Emitter2DFit{Float64} from (x, y, dataset) tuples,
# in microns. frame=1 throughout — the DBSCAN backend is frame-agnostic.
function _make_2d_smld(points::Vector{Tuple{Float64,Float64,Int}};
                      n_datasets::Int = maximum(p[3] for p in points))
    cam = IdealCamera(1:64, 1:64, 0.1)
    emitters = [Emitter2DFit{Float64}(
        x, y, 1000.0, 10.0, 0.01, 0.01, 50.0, 2.0;
        frame = 1, dataset = ds,
    ) for (x, y, ds) in points]
    BasicSMLD(emitters, cam, 1, n_datasets, Dict{String,Any}())
end

# Helper: generate `n` points in a Gaussian blob centered at (cx, cy) with σ in microns.
function _blob(rng, cx, cy, σ, n; dataset = 1)
    [(cx + σ * randn(rng), cy + σ * randn(rng), dataset) for _ in 1:n]
end

@testset "DBSCAN backend" begin

    @testset "config construction" begin
        cfg = DBSCANConfig(eps_nm = 50.0)
        @test cfg isa AbstractClusterConfig
        @test cfg.eps_nm == 50.0
        @test cfg.min_points == 5
        @test cfg.use_3d === false
        @test cfg.per_dataset === true
        @test cfg.remove_unclustered === false

        cfg2 = DBSCANConfig(eps_nm = 30.0, min_points = 3,
                            use_3d = true, per_dataset = false,
                            remove_unclustered = true)
        @test cfg2.min_points == 3
        @test cfg2.use_3d === true
        @test cfg2.per_dataset === false
        @test cfg2.remove_unclustered === true
    end

    @testset "three well-separated blobs + noise" begin
        rng = Xoshiro(20260417)
        σ = 0.010                 # 10 nm cluster scatter — tight
        n_per_blob = 40

        pts = Tuple{Float64,Float64,Int}[]
        # Centers 1 μm apart — far beyond eps=100 nm
        append!(pts, _blob(rng, 2.0, 2.0, σ, n_per_blob))
        append!(pts, _blob(rng, 4.0, 2.0, σ, n_per_blob))
        append!(pts, _blob(rng, 3.0, 4.0, σ, n_per_blob))
        # 30 noise points uniformly in [0, 6] × [0, 6], well outside any blob
        for _ in 1:30
            push!(pts, (6.0 * rand(rng), 6.0 * rand(rng), 1))
        end
        smld = _make_2d_smld(pts; n_datasets = 1)
        n_in = length(smld.emitters)

        cfg = DBSCANConfig(eps_nm = 100.0, min_points = 5, per_dataset = false)
        smld_out, info = cluster(smld, cfg)

        @test info isa ClusterInfo
        @test info.algorithm === :dbscan
        @test info.n_locs_in == n_in
        @test info.elapsed_s >= 0
        # With σ=10 nm and centers 1 μm apart, DBSCAN at eps=100 nm must find
        # exactly the three planted blobs.
        @test info.n_clusters == 3
        @test sum(info.cluster_sizes) == info.n_clustered
        @test info.n_clustered + info.n_noise == info.n_locs_in
        # All 120 blob points should cluster; almost all 30 noise points should
        # be flagged noise (generous bound accounts for occasional nearby draws).
        @test info.n_clustered >= 115
        @test info.n_noise >= 20
    end

    @testset "labels written to emitter.id + remove_unclustered + non-mutating" begin
        rng = Xoshiro(1)
        # Two tight blobs + 5 distant noise points.
        pts = Tuple{Float64,Float64,Int}[]
        append!(pts, _blob(rng, 1.0, 1.0, 0.005, 30))
        append!(pts, _blob(rng, 5.0, 5.0, 0.005, 30))
        for k in 1:5
            push!(pts, (10.0 + k, 10.0 + k, 1))
        end
        smld = _make_2d_smld(pts; n_datasets = 1)
        n_in = length(smld.emitters)

        # Default: remove_unclustered=false.
        cfg = DBSCANConfig(eps_nm = 100.0, min_points = 5, per_dataset = false)
        smld_keep, info_keep = cluster(smld, cfg)
        @test length(smld_keep.emitters) == n_in
        @test info_keep.n_clusters == 2
        # Input SMLD emitters are NOT mutated (non-mutating semantics).
        @test all(e -> e.id == 0, smld.emitters)
        # Output SMLD carries the labels.
        @test all(e -> e.id in 0:info_keep.n_clusters, smld_keep.emitters)
        @test any(e -> e.id == 0, smld_keep.emitters)
        # Each cluster's size should match what info reports.
        for k in 1:info_keep.n_clusters
            @test count(e -> e.id == k, smld_keep.emitters) == info_keep.cluster_sizes[k]
        end

        # remove_unclustered=true: same input SMLD, still not mutated.
        cfg_rm = DBSCANConfig(eps_nm = 100.0, min_points = 5,
                              per_dataset = false, remove_unclustered = true)
        smld_rm, info_rm = cluster(smld, cfg_rm)
        @test info_rm.n_clusters == 2
        @test length(smld_rm.emitters) == info_rm.n_clustered
        @test all(e -> e.id != 0, smld_rm.emitters)
        # Original smld still untouched.
        @test all(e -> e.id == 0, smld.emitters)
    end

    @testset "per_dataset label namespace is local" begin
        # Two datasets, each with two tight blobs; ids within each dataset should
        # be 1..2, so (dataset, id) is unique but id alone overlaps.
        rng = Xoshiro(2)
        pts = Tuple{Float64,Float64,Int}[]
        append!(pts, _blob(rng, 1.0, 1.0, 0.005, 20; dataset = 1))
        append!(pts, _blob(rng, 3.0, 3.0, 0.005, 20; dataset = 1))
        append!(pts, _blob(rng, 1.0, 1.0, 0.005, 20; dataset = 2))
        append!(pts, _blob(rng, 3.0, 3.0, 0.005, 20; dataset = 2))
        smld = _make_2d_smld(pts; n_datasets = 2)

        cfg = DBSCANConfig(eps_nm = 100.0, min_points = 5, per_dataset = true)
        smld_out, info = cluster(smld, cfg)
        @test info.n_clusters == 4
        @test info.n_clustered == 80

        # Within each dataset the ids should be {1, 2}.
        ids_ds1 = sort!(unique(e.id for e in smld_out.emitters if e.dataset == 1))
        ids_ds2 = sort!(unique(e.id for e in smld_out.emitters if e.dataset == 2))
        @test ids_ds1 == [1, 2]
        @test ids_ds2 == [1, 2]

        # Contrast: with per_dataset=false, all points cluster into 2 groups
        # (the two spatial centers are identical across datasets). Reuse the
        # same input smld — non-mutating semantics mean ids are still all 0.
        cfg_flat = DBSCANConfig(eps_nm = 100.0, min_points = 5, per_dataset = false)
        _, info_flat = cluster(smld, cfg_flat)
        @test info_flat.n_clusters == 2
        @test info_flat.n_clustered == 80
    end

    @testset "argument validation" begin
        smld = _make_2d_smld([(0.0, 0.0, 1)])
        @test_throws ArgumentError cluster(smld, DBSCANConfig(eps_nm = 0.0))
        @test_throws ArgumentError cluster(smld, DBSCANConfig(eps_nm = 50.0, min_points = 0))
    end

    @testset "use_3d on 2D data errors" begin
        smld = _make_2d_smld([(0.0, 0.0, 1), (1.0, 1.0, 1)])
        cfg = DBSCANConfig(eps_nm = 100.0, use_3d = true)
        @test_throws ErrorException cluster(smld, cfg)
    end

    @testset "3D clustering path" begin
        cam = IdealCamera(1:64, 1:64, 0.1)
        rng = Xoshiro(3)
        emitters = SMLMData.Emitter3DFit{Float64}[]
        # Two tight 3D blobs separated in z only (centers (2,2,0) and (2,2,1) μm).
        for _ in 1:25
            push!(emitters, Emitter3DFit{Float64}(
                2.0 + 0.005 * randn(rng), 2.0 + 0.005 * randn(rng), 0.0 + 0.005 * randn(rng),
                1000.0, 10.0, 0.01, 0.01, 0.02, 50.0, 2.0; frame=1, dataset=1))
        end
        for _ in 1:25
            push!(emitters, Emitter3DFit{Float64}(
                2.0 + 0.005 * randn(rng), 2.0 + 0.005 * randn(rng), 1.0 + 0.005 * randn(rng),
                1000.0, 10.0, 0.01, 0.01, 0.02, 50.0, 2.0; frame=1, dataset=1))
        end
        smld3 = BasicSMLD(emitters, cam, 1, 1, Dict{String,Any}())

        cfg = DBSCANConfig(eps_nm = 100.0, min_points = 5, use_3d = true, per_dataset = false)
        _, info = cluster(smld3, cfg)
        @test info.n_clusters == 2
        @test info.n_clustered == 50
    end

    @testset "empty SMLD" begin
        smld = _make_2d_smld(Tuple{Float64,Float64,Int}[]; n_datasets = 1)
        cfg = DBSCANConfig(eps_nm = 50.0, per_dataset = false)
        smld_out, info = cluster(smld, cfg)
        @test info.n_locs_in == 0
        @test info.n_clusters == 0
        @test info.n_clustered == 0
        @test info.n_noise == 0
        @test isempty(info.cluster_sizes)
        @test isempty(smld_out.emitters)
    end

end
