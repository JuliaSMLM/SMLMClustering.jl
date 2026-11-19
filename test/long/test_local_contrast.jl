using SMLMClustering
using SMLMData
using Test
using Random
using Statistics
include(joinpath(@__DIR__, "utils", "helpers.jl"))

@testset "LocalContrastFeature (thorough)" begin

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

    @testset "global density gradient: contrast cancels gradient" begin
        rng = Xoshiro(20260502)
        # Linearly-increasing density along x: a small region in the
        # low-density (left) side that is locally elevated should get higher
        # contrast than the globally-denser right side, even though the right
        # side has higher *absolute* density.
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

end
