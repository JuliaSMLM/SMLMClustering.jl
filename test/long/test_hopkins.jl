using SMLMClustering
using SMLMData
using Test
using Random
include(joinpath(@__DIR__, "utils", "helpers.jl"))

@testset "Hopkins backend (thorough)" begin

    @testset "uniform random data → H ≈ 0.5" begin
        # Uniform-random 2D over a 5×5 μm box, single dataset; with enough
        # samples and repeats the expected H is 0.5. Allow ±0.1 slack.
        rng = Xoshiro(20260427)
        n = 2000
        pts = Tuple{Float64,Float64,Int}[]
        for _ in 1:n
            push!(pts, (5.0 * rand(rng), 5.0 * rand(rng), 1))
        end
        smld = _make_2d_smld(pts; n_datasets = 1)

        cfg = HopkinsConfig(n_samples = 100, random_repeats = 10, seed = 1,
                            per_dataset = false)
        smld_out, info = cluster_statistics(smld, cfg)
        @test smld_out === smld    # passthrough — same reference, no copy
        @test info isa ClusterStatisticsInfo
        @test info.algorithm === :hopkins
        @test info.statistic_name === :hopkins
        @test info.n_locs_in == n
        @test info.elapsed_s >= 0
        @test 0.4 <= info.statistic <= 0.6
    end

    @testset "tight gaussian blob → H > 0.85" begin
        rng = Xoshiro(20260427)
        # Single very tight blob of 500 points (σ = 5 nm) — strong clustering.
        pts = Tuple{Float64,Float64,Int}[]
        for _ in 1:500
            push!(pts, (2.5 + 0.005 * randn(rng), 2.5 + 0.005 * randn(rng), 1))
        end
        # A handful of scattered points to give the bbox some extent
        # (otherwise the bbox collapses to ~σ and the test loses meaning).
        for k in 1:8
            push!(pts, (rand(rng) * 5.0, rand(rng) * 5.0, 1))
        end
        smld = _make_2d_smld(pts; n_datasets = 1)

        cfg = HopkinsConfig(n_samples = 50, random_repeats = 10, seed = 2,
                            per_dataset = false)
        _, info = cluster_statistics(smld, cfg)
        @test info.statistic > 0.85
    end

    @testset "per_dataset produces per-dataset vector" begin
        rng = Xoshiro(20260427)
        # Two datasets with very different clustering tendency:
        # ds 1: tight blob (H high)  /  ds 2: uniform (H ≈ 0.5)
        pts = Tuple{Float64,Float64,Int}[]
        for _ in 1:500
            push!(pts, (2.5 + 0.005 * randn(rng), 2.5 + 0.005 * randn(rng), 1))
        end
        # Scatter a few points in ds 1 to extend the bbox
        for _ in 1:8
            push!(pts, (5.0 * rand(rng), 5.0 * rand(rng), 1))
        end
        for _ in 1:1000
            push!(pts, (5.0 * rand(rng), 5.0 * rand(rng), 2))
        end
        smld = _make_2d_smld(pts; n_datasets = 2)

        cfg = HopkinsConfig(n_samples = 50, random_repeats = 5, seed = 3,
                            per_dataset = true)
        _, info = cluster_statistics(smld, cfg)
        per_ds = info.extras[:hopkins_per_dataset]
        @test per_ds isa Vector{Float64}
        @test length(per_ds) == 2
        @test per_ds[1] > 0.8           # tight blob
        @test 0.4 <= per_ds[2] <= 0.6   # uniform
        # Reported statistic is the mean across datasets.
        @test info.statistic ≈ (per_ds[1] + per_ds[2]) / 2
    end

    @testset "seed reproducibility" begin
        rng = Xoshiro(20260427)
        pts = Tuple{Float64,Float64,Int}[]
        for _ in 1:1000
            push!(pts, (5.0 * rand(rng), 5.0 * rand(rng), 1))
        end
        smld = _make_2d_smld(pts; n_datasets = 1)

        cfg = HopkinsConfig(n_samples = 50, random_repeats = 3, seed = 99,
                            per_dataset = false)
        _, info_a = cluster_statistics(smld, cfg)
        _, info_b = cluster_statistics(smld, cfg)
        @test info_a.statistic == info_b.statistic

        # Different seed → almost certainly different result on uniform data.
        cfg2 = HopkinsConfig(n_samples = 50, random_repeats = 3, seed = 100,
                             per_dataset = false)
        _, info_c = cluster_statistics(smld, cfg2)
        @test info_a.statistic != info_c.statistic
    end

    @testset "3D data path" begin
        cam = IdealCamera(1:64, 1:64, 0.1)
        rng = Xoshiro(20260427)
        emitters = SMLMData.Emitter3DFit{Float64}[]
        # 3D uniform random in a 5×5×5 μm box.
        for _ in 1:1000
            push!(emitters, Emitter3DFit{Float64}(
                5.0 * rand(rng), 5.0 * rand(rng), 5.0 * rand(rng),
                1000.0, 10.0, 0.01, 0.01, 0.02, 50.0, 2.0; frame = 1, dataset = 1))
        end
        smld3 = BasicSMLD(emitters, cam, 1, 1, Dict{String,Any}())

        cfg = HopkinsConfig(n_samples = 80, random_repeats = 8, seed = 5,
                            use_3d = true, per_dataset = false)
        _, info = cluster_statistics(smld3, cfg)
        @test 0.4 <= info.statistic <= 0.6
    end

    @testset "empty SMLD returns NaN, not an error" begin
        smld = _make_2d_smld(Tuple{Float64,Float64,Int}[]; n_datasets = 1)
        cfg = HopkinsConfig(n_samples = 10, per_dataset = false)
        _, info = cluster_statistics(smld, cfg)
        @test info.n_locs_in == 0
        @test isnan(info.statistic)
    end

    @testset "n_samples > n_points returns NaN per group" begin
        # Tiny SMLD: only 5 points, n_samples=10 → too many samples requested.
        smld = _make_2d_smld([(0.0, 0.0, 1), (1.0, 0.5, 1), (0.3, 1.2, 1),
                              (1.4, 1.4, 1), (0.7, 0.2, 1)])
        cfg = HopkinsConfig(n_samples = 10, per_dataset = false)
        _, info = cluster_statistics(smld, cfg)
        @test isnan(info.statistic)
    end

    @testset "region restricts reference points to the polygon" begin
        # Data uniform INSIDE a non-convex U-shape. bbox-sampled references land in
        # the empty central notch → H falsely high; restricting references to the
        # polygon recovers the random null (~0.5). Regression for the window bug.
        poly = [(0.0, 0.0), (5.0, 0.0), (5.0, 5.0), (4.0, 5.0),
                (4.0, 1.0), (1.0, 1.0), (1.0, 5.0), (0.0, 5.0)]
        rng = Xoshiro(42)
        pts = Tuple{Float64,Float64,Int}[]
        while length(pts) < 2000
            x = 5rand(rng); y = 5rand(rng)
            SMLMClustering._point_in_polygon(x, y, poly) && push!(pts, (x, y, 1))
        end
        smld = _make_2d_smld(pts; n_datasets = 1)

        base(; kw...) = HopkinsConfig(; n_samples = 150, random_repeats = 15,
                                      seed = 1, per_dataset = false, kw...)
        H_bbox   = cluster_statistics(smld, base())[2].statistic
        H_region = cluster_statistics(smld, base(region = poly))[2].statistic
        @test H_bbox > 0.6                 # falsely clustered without a region
        @test 0.4 <= H_region <= 0.6       # random null once references are confined
        @test H_bbox > H_region + 0.15     # relational, robust to constants

        # metadata channel: the same polygon supplied via smld.metadata
        smld_meta = SMLMData.BasicSMLD(smld.emitters, smld.camera, smld.n_frames,
                                       smld.n_datasets,
                                       Dict{String,Any}("edge_outer_polygon" => poly))
        H_meta = cluster_statistics(smld_meta, base(region = :metadata))[2].statistic
        @test 0.4 <= H_meta <= 0.6
    end

end
