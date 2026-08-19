using SMLMClustering
using SMLMData
using SMLMSim
using Random
using Test

@testset "boundary_statistics backend" begin

    if SMLM_TEST_FULL
    @testset "boundary_statistics — SMLMSim integration: simulate → cluster → boundary_clusters → boundary_statistics" begin
        # Simulate a realistic small set of localizations from a ground
        # truth of several Nmer clusters (SMLMSim), instead of loading a
        # real dataset.
        Random.seed!(42)   # reproducible run
        cam         = IdealCamera(1:40, 1:40, 0.1)   # 4 μm × 4 μm field of view
        sim_pattern = Nmer2D(n = 6, d = 0.1)         # 6-molecule cluster, 100 nm diameter
        sim_params  = StaticSMLMConfig(
            density    = 0.5,   # patterns/μm² → ~8 clusters over the 16 μm² field
            σ_psf      = 0.02,  # 20 nm localization precision
            minphotons = 50,
            ndatasets  = 1,
            nframes    = 300,
            framerate  = 50.0,
            ndims      = 2,
        )
        # Higher k_on than the package default so blinking yields a
        # reasonable number of localizations over just 300 frames (default
        # k_on=1e-2 gives an equilibrium on-fraction of ~0.02%, i.e. only a
        # handful of localizations total; k_on=1.0 gives ~2%, i.e. hundreds
        # — a realistic "small set").
        sim_molecule  = GenericFluor(photons = 1e4, k_off = 50.0, k_on = 1.0)
        smld, siminfo = simulate(sim_params; pattern = sim_pattern, camera = cam,
                                  molecule = sim_molecule)

        per_dataset = false

        cfg_dbscan = DBSCANConfig(
            eps_nm             = 150.0, # neighborhood radius in nm (required)
            min_points         = 5,     # core-point threshold / min cluster size
            use_3d             = false, # include z-coordinate
            per_dataset        = per_dataset,
            remove_unclustered = false,
        )
        (smld_out, info_dbscan) = cluster(smld, cfg_dbscan)

        cfg_bnd = BoundaryClustersConfig(
            shrink      = 0.05,  # 0 = convex hull, 1 = tightest alpha-shape
            per_dataset = per_dataset,
        )
        (_, binfo) = boundary_clusters(smld_out, cfg_bnd)

        (_, sinfo) = boundary_statistics(smld_out, binfo; per_dataset = per_dataset,
                                          algorithm = info_dbscan.algorithm)

        @test siminfo.n_patterns >= 3          # "several" ground-truth clusters
        @test siminfo.n_localizations >= 20    # realistic small-but-nontrivial set
        @test info_dbscan.n_clusters >= 1
        @test sinfo isa ClusterBoundaryStatisticsInfo
        @test sinfo.n_clusters == binfo.n_clusters == info_dbscan.n_clusters
        @test sinfo.cluster_keys == binfo.cluster_keys
        @test sinfo.n_locs_in == length(smld.emitters)
        @test sinfo.n_clustered + sinfo.n_noise == sinfo.n_locs_in
        @test sinfo.algorithm == :dbscan
        @test sinfo.shrink == cfg_bnd.shrink
        @test all(a -> a > 0, sinfo.area)
        @test all(p -> p > 0, sinfo.perimeter)
        @test all(sinfo.n_pts_per_cluster .>= cfg_dbscan.min_points)
        @test sum(sinfo.n_pts_per_cluster) == sinfo.n_clustered
        @test all(c -> all(isfinite, c), sinfo.center)
        # compactness = 4π·area/perimeter² is bounded in (0,1] by the
        # isoperimetric inequality for any simple closed polygon — always safe.
        @test all(x -> 0.0 < x <= 1.0 + 1e-9, sinfo.compactness)
        # circularity/solidity/convexity depend on the shrink=1 "tightest"
        # alpha-shape boundary, which can collapse to a small sub-loop for
        # dense clusters (same documented quirk as test_boundary.jl's
        # convexity note) — assert finiteness/sign, not a (0,1] bound.
        @test all(x -> isfinite(x) && x >= 0, sinfo.circularity)
        @test all(x -> isfinite(x) && x >= 0, sinfo.solidity)
        @test all(x -> isfinite(x) && x > 0, sinfo.convexity)
        @test length(sinfo.nn_within_clusters) == sinfo.n_clusters
        @test all(j -> length(sinfo.nn_within_clusters[j]) == sinfo.n_pts_per_cluster[j],
                  1:sinfo.n_clusters)
        @test isfinite(sinfo.min_c2c_dist)   # multiple clusters ⇒ well-defined
        @test isfinite(sinfo.min_e2e_dist)
        @test sinfo.elapsed_s >= 0
    end
    end  # SMLM_TEST_FULL

end
