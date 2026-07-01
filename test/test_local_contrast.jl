using SMLMClustering
using SMLMData
using Test
using Random
using Statistics

# Reuses `_make_2d_smld` and `_blob` from test_dbscan.jl (top-level).

@testset "LocalContrastFeature" begin

    @testset "config construction + validation" begin
        cfg = LocalContrastFeature()
        @test cfg isa AbstractStatisticsConfig
        @test cfg.density_k == 200
        @test cfg.background_k == 2000
        @test cfg.use_3d === false
        @test cfg.per_dataset === false

        cfg2 = LocalContrastFeature(density_k = 20, background_k = 200,
                                    per_dataset = true)
        @test cfg2.density_k == 20
        @test cfg2.background_k == 200
        @test cfg2.per_dataset === true

        # background_k must be > density_k
        rng = Xoshiro(20260502)
        pts = [(0.5 * randn(rng), 0.5 * randn(rng), 1) for _ in 1:50]
        smld = _make_2d_smld(pts; n_datasets = 1)
        @test_throws ArgumentError cluster_statistics(
            smld, LocalContrastFeature(density_k = 50, background_k = 50))
        @test_throws ArgumentError cluster_statistics(
            smld, LocalContrastFeature(density_k = 100, background_k = 50))
        @test_throws ArgumentError cluster_statistics(
            smld, LocalContrastFeature(density_k = 0, background_k = 10))
    end

    @testset "shape and pass-through semantics" begin
        rng = Xoshiro(20260502)
        pts = [(0.5 * randn(rng), 0.5 * randn(rng), 1) for _ in 1:200]
        smld = _make_2d_smld(pts; n_datasets = 1)

        cfg = LocalContrastFeature(density_k = 10, background_k = 50)
        smld_out, info = cluster_statistics(smld, cfg)

        @test smld_out === smld    # pass-through (V10)
        @test info isa ClusterStatisticsInfo
        @test info.algorithm === :local_contrast
        @test info.statistic_name === :median_local_contrast
        @test info.n_locs_in == 200
        @test info.elapsed_s >= 0

        c = info.extras[:contrast_per_emitter]
        f = info.extras[:log_density_per_emitter]
        @test c isa Vector{Float64}
        @test f isa Vector{Float64}
        @test length(c) == 200
        @test length(f) == 200
    end

    @testset "small group → NaN feature when n ≤ density_k" begin
        rng = Xoshiro(20260502)
        # 8 points, density_k=10: cannot compute → NaN entries.
        pts = [(0.1 * randn(rng), 0.1 * randn(rng), 1) for _ in 1:8]
        smld = _make_2d_smld(pts; n_datasets = 1)
        cfg = LocalContrastFeature(density_k = 10, background_k = 20)
        _, info = cluster_statistics(smld, cfg)
        @test all(isnan, info.extras[:contrast_per_emitter])
        @test all(isnan, info.extras[:log_density_per_emitter])
        @test isnan(info.statistic)
    end

    if SMLM_TEST_FULL
    @testset "tight blob in sparse field → contrast positive in blob, ≈0 elsewhere" begin
        rng = Xoshiro(20260502)
        # 1500 uniform background + 200-point tight blob at (3, 3).
        # Blob points should have contrast >> 0; background contrast should
        # cluster around 0.
        pts = [(5.0 * rand(rng), 5.0 * rand(rng), 1) for _ in 1:1500]
        n_blob = 200
        blob_start = length(pts) + 1
        append!(pts, _blob(rng, 3.0, 3.0, 0.05, n_blob))
        smld = _make_2d_smld(pts; n_datasets = 1)

        cfg = LocalContrastFeature(density_k = 20, background_k = 400)
        _, info = cluster_statistics(smld, cfg)
        c = info.extras[:contrast_per_emitter]

        blob_idxs = blob_start:length(pts)
        bg_idxs = 1:(blob_start - 1)

        blob_mean = mean(filter(isfinite, c[blob_idxs]))
        bg_median = median(filter(isfinite, c[bg_idxs]))

        @test blob_mean > 1.0          # well above ≥0.25 typical seed cutoff
        @test abs(bg_median) < 0.3     # background contrast clusters around 0
        @test blob_mean > bg_median + 1.0
    end
    end  # SMLM_TEST_FULL

    if SMLM_TEST_FULL
    @testset "global density gradient: contrast cancels gradient" begin
        rng = Xoshiro(20260502)
        # Linearly-increasing density along x: a small region in the
        # low-density (left) side that is locally elevated should get higher
        # contrast than the globally-denser right side, even though the right
        # side has higher *absolute* density. This is the RGY left/right fix.
        pts = Tuple{Float64, Float64, Int}[]
        # Background: density rises with x (sample x from triangular distribution).
        for _ in 1:3000
            x = 5.0 * sqrt(rand(rng))   # CDF gives linear ramp 0..5
            y = 5.0 * rand(rng)
            push!(pts, (x, y, 1))
        end
        # Local elevation on the LEFT side: small dense blob.
        local_blob_start = length(pts) + 1
        append!(pts, _blob(rng, 1.0, 2.5, 0.05, 150))
        right_idxs = [i for (i, p) in pairs(pts) if i < local_blob_start && p[1] > 4.0]

        smld = _make_2d_smld(pts; n_datasets = 1)
        cfg = LocalContrastFeature(density_k = 20, background_k = 400)
        _, info = cluster_statistics(smld, cfg)
        c = info.extras[:contrast_per_emitter]
        f = info.extras[:log_density_per_emitter]

        # Absolute density on the right side is HIGHER than at the local blob
        # on the left, but local contrast at the blob is HIGHER than on the
        # globally-denser right side.
        right_logden = mean(filter(isfinite, f[right_idxs]))
        blob_logden = mean(filter(isfinite, f[local_blob_start:end]))
        right_contrast = mean(filter(isfinite, c[right_idxs]))
        blob_contrast = mean(filter(isfinite, c[local_blob_start:end]))

        @test blob_logden > right_logden             # blob has higher absolute density too
        @test blob_contrast > right_contrast + 0.5   # but contrast separates much more
        @test right_contrast < 0.3                   # gradient is cancelled
    end
    end  # SMLM_TEST_FULL

    @testset "per_dataset isolates feature computation" begin
        rng = Xoshiro(20260502)
        pts = Tuple{Float64, Float64, Int}[]
        # Two datasets at the SAME location with different point counts.
        # Pooled: bg=100 spans the union → sees ds1 + ds2 jointly.
        # Per-dataset: bg=100 is per group → only sees own dataset.
        # The densities differ in this overlap regime.
        for _ in 1:200
            push!(pts, (0.5 * rand(rng), 0.5 * rand(rng), 1))
        end
        for _ in 1:200
            push!(pts, (0.5 * rand(rng), 0.5 * rand(rng), 2))
        end
        smld = _make_2d_smld(pts; n_datasets = 2)

        _, info_pooled = cluster_statistics(
            smld, LocalContrastFeature(density_k = 20, background_k = 100,
                                       per_dataset = false))
        _, info_per = cluster_statistics(
            smld, LocalContrastFeature(density_k = 20, background_k = 100,
                                       per_dataset = true))
        # Log-density at fine scale: pooled sees ~2x the density (both datasets
        # in the same area), so pooled log-densities are ≈ +log(2) higher than
        # per-dataset for the same point.
        f_pool = info_pooled.extras[:log_density_per_emitter]
        f_per = info_per.extras[:log_density_per_emitter]
        @test f_pool != f_per
        # Median pooled log-density should be ~log(2) ≈ 0.69 above per-dataset.
        @test median(filter(isfinite, f_pool)) >
              median(filter(isfinite, f_per)) + 0.4
    end

end  # @testset LocalContrastFeature
