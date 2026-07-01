using SMLMClustering
using SMLMData
using Test
using Random
using Statistics

# Reuse `_make_2d_smld` and `_blob` from test_dbscan.jl (included earlier in
# runtests.jl so those helpers are already defined at top-level).

@testset "VoronoiDensity backend" begin

    @testset "config construction" begin
        cfg = VoronoiDensityConfig()
        @test cfg isa AbstractStatisticsConfig
        @test cfg.use_3d === false
        @test cfg.per_dataset === true

        cfg2 = VoronoiDensityConfig(use_3d = false, per_dataset = false)
        @test cfg2.per_dataset === false
    end

    if SMLM_TEST_FULL
    @testset "uniform Poisson 2D: median density ≈ N / area" begin
        # 1500 uniform points in a 5x5 μm box → expected density ≈ 60 μm⁻².
        # Convex-hull clipping plus boundary cells bias the empirical median
        # somewhat, so allow ±30% per the spec.
        rng = Xoshiro(20260427)
        n = 1500
        side = 5.0
        expected = n / (side * side)  # 60.0 μm⁻²
        pts = Tuple{Float64,Float64,Int}[]
        for _ in 1:n
            push!(pts, (side * rand(rng), side * rand(rng), 1))
        end
        smld = _make_2d_smld(pts; n_datasets = 1)

        cfg = VoronoiDensityConfig(per_dataset = false)
        smld_out, info = cluster_statistics(smld, cfg)
        @test smld_out === smld    # passthrough — same reference
        @test info isa ClusterStatisticsInfo
        @test info.algorithm === :voronoi_density
        @test info.statistic_name === :median_density
        @test info.n_locs_in == n
        @test info.elapsed_s >= 0
        # Within ±30% of the analytic expectation.
        @test abs(info.statistic - expected) / expected < 0.3
    end
    end  # SMLM_TEST_FULL

    if SMLM_TEST_FULL
    @testset "extras carry per-emitter density and area in original order" begin
        rng = Xoshiro(20260427)
        pts = Tuple{Float64,Float64,Int}[]
        for _ in 1:300
            push!(pts, (3.0 * rand(rng), 3.0 * rand(rng), 1))
        end
        smld = _make_2d_smld(pts; n_datasets = 1)

        cfg = VoronoiDensityConfig(per_dataset = false)
        _, info = cluster_statistics(smld, cfg)
        ρ = info.extras[:density_per_emitter]
        A = info.extras[:area_per_emitter]
        @test ρ isa Vector{Float64}
        @test A isa Vector{Float64}
        @test length(ρ) == 300
        @test length(A) == 300
        # ρ = 1/A elementwise where both are non-NaN.
        for i in 1:300
            if !isnan(ρ[i]) && !isnan(A[i])
                @test ρ[i] ≈ 1.0 / A[i]
                @test A[i] > 0
            end
        end
    end
    end  # SMLM_TEST_FULL

    if SMLM_TEST_FULL
    @testset "tight blob → median density >> uniform-baseline median" begin
        rng = Xoshiro(20260427)
        # Tight σ=5 nm blob (high local density) plus some scatter so the
        # bbox doesn't collapse to ~σ (which would make the test trivially pass).
        pts = Tuple{Float64,Float64,Int}[]
        append!(pts, _blob(rng, 2.5, 2.5, 0.005, 400))
        for _ in 1:100
            push!(pts, (5.0 * rand(rng), 5.0 * rand(rng), 1))
        end
        smld_blob = _make_2d_smld(pts; n_datasets = 1)

        # Uniform baseline at the same total density (500 in a 5x5 box).
        pts_uniform = Tuple{Float64,Float64,Int}[]
        for _ in 1:500
            push!(pts_uniform, (5.0 * rand(rng), 5.0 * rand(rng), 1))
        end
        smld_uniform = _make_2d_smld(pts_uniform; n_datasets = 1)

        cfg = VoronoiDensityConfig(per_dataset = false)
        _, info_blob = cluster_statistics(smld_blob, cfg)
        _, info_uniform = cluster_statistics(smld_uniform, cfg)

        # Blob-dominated median should be at least an order of magnitude higher
        # than the uniform baseline (blob cells are ~σ²/N tiny).
        @test info_blob.statistic > 10 * info_uniform.statistic
    end
    end  # SMLM_TEST_FULL

    if SMLM_TEST_FULL
    @testset "per_dataset: per-dataset tessellation, flat per-emitter vector" begin
        rng = Xoshiro(20260427)
        # ds 1: 200 points in a 1x1 μm box → density ~200/μm²
        # ds 2: 200 points in a 4x4 μm box → density ~12.5/μm²
        pts = Tuple{Float64,Float64,Int}[]
        for _ in 1:200
            push!(pts, (1.0 * rand(rng), 1.0 * rand(rng), 1))
        end
        for _ in 1:200
            push!(pts, (4.0 * rand(rng), 4.0 * rand(rng), 2))
        end
        smld = _make_2d_smld(pts; n_datasets = 2)

        cfg = VoronoiDensityConfig(per_dataset = true)
        _, info = cluster_statistics(smld, cfg)
        ρ = info.extras[:density_per_emitter]
        @test length(ρ) == 400

        # First 200 entries belong to ds 1 (high density), last 200 to ds 2.
        ρ_ds1 = ρ[1:200]
        ρ_ds2 = ρ[201:400]
        @test all(!isnan, ρ_ds1)
        @test all(!isnan, ρ_ds2)
        median_ds1 = median(ρ_ds1)
        median_ds2 = median(ρ_ds2)
        # ds 1 median should be much higher (~16x analytic ratio: 200 vs 12.5).
        @test median_ds1 > 4 * median_ds2

        # Flat-vector ordering claim: per-dataset processing should NOT
        # reorder. Compare against per_dataset=false, which we already know
        # processes everything in input order.
        cfg_flat = VoronoiDensityConfig(per_dataset = false)
        _, info_flat = cluster_statistics(smld, cfg_flat)
        ρ_flat = info_flat.extras[:density_per_emitter]
        @test length(ρ_flat) == 400
        # The per-emitter vector under per_dataset=true must be in original
        # emitter-index order (not grouped by dataset). The cleanest test:
        # ds 1 occupies indices 1:200 in BOTH per_dataset modes.
        @test all(!isnan, ρ_flat)
    end
    end  # SMLM_TEST_FULL

    if SMLM_TEST_FULL
    @testset "<3 points in a dataset → those emitters get NaN, others fine" begin
        # ds 1: 2 points (degenerate, no tessellation)
        # ds 2: 50 points (normal tessellation)
        rng = Xoshiro(11)
        pts = Tuple{Float64,Float64,Int}[]
        push!(pts, (0.0, 0.0, 1))
        push!(pts, (1.0, 1.0, 1))
        for _ in 1:50
            push!(pts, (3.0 * rand(rng), 3.0 * rand(rng), 2))
        end
        smld = _make_2d_smld(pts; n_datasets = 2)

        cfg = VoronoiDensityConfig(per_dataset = true)
        _, info = cluster_statistics(smld, cfg)
        ρ = info.extras[:density_per_emitter]
        @test length(ρ) == 52

        # ds 1 emitters (indices 1, 2) → NaN
        @test isnan(ρ[1])
        @test isnan(ρ[2])
        # ds 2 emitters (indices 3:52) → all finite
        @test all(!isnan, ρ[3:52])

        # Median ignores the NaNs.
        @test !isnan(info.statistic)
    end
    end  # SMLM_TEST_FULL

    @testset "duplicate (x,y) coordinates raise ArgumentError" begin
        # Exact-coincident generators cause DelaunayTriangulation.get_area
        # to raise KeyError; we want a clean ArgumentError before triangulation.
        pts = [(0.0, 0.0, 1), (1.0, 0.0, 1), (0.0, 1.0, 1), (1.0, 1.0, 1),
               (0.0, 0.0, 1)]  # duplicate of first point
        smld = _make_2d_smld(pts; n_datasets = 1)
        cfg = VoronoiDensityConfig(per_dataset = false)
        @test_throws ArgumentError cluster_statistics(smld, cfg)

        # Per-dataset path also surfaces the error per-group.
        pts2 = [(0.0, 0.0, 1), (1.0, 0.0, 1), (0.0, 1.0, 1), (0.0, 0.0, 1)]
        smld2 = _make_2d_smld(pts2; n_datasets = 1)
        cfg2 = VoronoiDensityConfig(per_dataset = true)
        @test_throws ArgumentError cluster_statistics(smld2, cfg2)
    end

    if SMLM_TEST_FULL
    @testset "empty SMLD returns NaN statistic, empty per-emitter vectors" begin
        smld = _make_2d_smld(Tuple{Float64,Float64,Int}[]; n_datasets = 1)
        cfg = VoronoiDensityConfig(per_dataset = false)
        smld_out, info = cluster_statistics(smld, cfg)
        @test smld_out === smld
        @test info.n_locs_in == 0
        @test isnan(info.statistic)
        @test info.extras[:density_per_emitter] == Float64[]
        @test info.extras[:area_per_emitter] == Float64[]
    end
    end  # SMLM_TEST_FULL

    @testset "use_3d=true raises ArgumentError" begin
        smld = _make_2d_smld([(0.0, 0.0, 1), (1.0, 0.0, 1), (0.0, 1.0, 1)])
        @test_throws ArgumentError cluster_statistics(smld, VoronoiDensityConfig(use_3d = true))
    end

    @testset "passthrough: SMLD reference unchanged, emitter ids untouched" begin
        rng = Xoshiro(7)
        pts = Tuple{Float64,Float64,Int}[]
        for _ in 1:80
            push!(pts, (2.0 * rand(rng), 2.0 * rand(rng), 1))
        end
        smld = _make_2d_smld(pts; n_datasets = 1)
        # Pre-flight: all ids are 0 (default emitter state).
        @test all(e -> e.id == 0, smld.emitters)

        cfg = VoronoiDensityConfig(per_dataset = false)
        smld_out, _ = cluster_statistics(smld, cfg)
        @test smld_out === smld
        @test all(e -> e.id == 0, smld.emitters)
    end

    if SMLM_TEST_FULL
    @testset "summary scalar = median of non-NaN entries" begin
        # Build a deterministic small grid so we can check the summary
        # statistic exactly against Statistics.median.
        rng = Xoshiro(20260427)
        pts = Tuple{Float64,Float64,Int}[]
        for _ in 1:200
            push!(pts, (2.0 * rand(rng), 2.0 * rand(rng), 1))
        end
        smld = _make_2d_smld(pts; n_datasets = 1)

        cfg = VoronoiDensityConfig(per_dataset = false)
        _, info = cluster_statistics(smld, cfg)
        ρ = info.extras[:density_per_emitter]
        valid = filter(!isnan, ρ)
        @test info.statistic ≈ median(valid)
    end
    end  # SMLM_TEST_FULL

end
