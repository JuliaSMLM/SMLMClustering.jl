using SMLMClustering
using SMLMData
using Test
using Random
using Statistics
include(joinpath(@__DIR__, "long", "utils", "helpers.jl"))

@testset "MRFDensityCluster backend" begin

    @testset "config construction" begin
        cfg = MRFDensityClusterConfig()
        @test cfg isa AbstractClusterConfig
        @test cfg.n_regimes == 2
        @test cfg.regime_thresholds === nothing
        @test cfg.regime_gaussians === nothing
        @test cfg.density_estimator === :voronoi
        @test cfg.density_k == 20
        @test cfg.smoothness_lambda === nothing
        @test cfg.graph_kind === :delaunay
        @test cfg.graph_k == 8
        @test cfg.inference === :icm
        @test cfg.icm_iters == 50
        @test cfg.min_points == 5
        @test cfg.use_3d === false
        @test cfg.per_dataset === true
        @test cfg.remove_unclustered === false

        cfg2 = MRFDensityClusterConfig(n_regimes = 3, graph_kind = :knn,
                                       graph_k = 12, min_points = 10,
                                       smoothness_lambda = 0.5,
                                       per_dataset = false,
                                       remove_unclustered = true)
        @test cfg2.n_regimes == 3
        @test cfg2.graph_kind === :knn
        @test cfg2.graph_k == 12
        @test cfg2.min_points == 10
        @test cfg2.smoothness_lambda == 0.5
        @test cfg2.per_dataset === false
        @test cfg2.remove_unclustered === true

        cfg3 = MRFDensityClusterConfig(n_regimes = 3,
                                       regime_thresholds = [3.0, 5.0])
        @test cfg3.regime_thresholds == [3.0, 5.0]

        cfg4 = MRFDensityClusterConfig(density_estimator = :knn, density_k = 30)
        @test cfg4.density_estimator === :knn
        @test cfg4.density_k == 30

        gaussians = (means = [2.0, 5.0], vars = [0.5, 0.5], weights = [0.7, 0.3])
        cfg5 = MRFDensityClusterConfig(regime_gaussians = gaussians)
        @test cfg5.regime_gaussians == gaussians
    end

    @testset "argument validation" begin
        # 3-point dummy SMLD so we get past the n<3 short-circuit and reach validation.
        smld_dummy = _make_2d_smld([(0.0, 0.0, 1), (1.0, 0.0, 1), (0.0, 1.0, 1)])

        @test_throws ArgumentError cluster(smld_dummy, MRFDensityClusterConfig(n_regimes = 1))
        @test_throws ArgumentError cluster(smld_dummy, MRFDensityClusterConfig(graph_kind = :bogus))
        @test_throws ArgumentError cluster(smld_dummy, MRFDensityClusterConfig(inference = :graph_cuts))
        @test_throws ArgumentError cluster(smld_dummy, MRFDensityClusterConfig(graph_k = 0))
        @test_throws ArgumentError cluster(smld_dummy, MRFDensityClusterConfig(density_estimator = :bogus))
        @test_throws ArgumentError cluster(smld_dummy, MRFDensityClusterConfig(density_k = 0))
        @test_throws ArgumentError cluster(smld_dummy, MRFDensityClusterConfig(icm_iters = 0))
        @test_throws ArgumentError cluster(smld_dummy, MRFDensityClusterConfig(min_points = 0))
        @test_throws ArgumentError cluster(smld_dummy, MRFDensityClusterConfig(smoothness_lambda = -1.0))
        # Threshold length mismatch (n_regimes=2 expects length 1, given 2).
        @test_throws ArgumentError cluster(smld_dummy,
            MRFDensityClusterConfig(n_regimes = 2, regime_thresholds = [1.0, 2.0]))
        # Unsorted thresholds.
        @test_throws ArgumentError cluster(smld_dummy,
            MRFDensityClusterConfig(n_regimes = 3, regime_thresholds = [5.0, 3.0]))
        # Calibration override modes are mutually exclusive.
        @test_throws ArgumentError cluster(smld_dummy,
            MRFDensityClusterConfig(regime_thresholds = [1.0],
                                    regime_gaussians = (means = [0.0, 1.0],
                                                        vars = [1.0, 1.0],
                                                        weights = [0.5, 0.5])))
        @test_throws ArgumentError cluster(smld_dummy,
            MRFDensityClusterConfig(regime_gaussians = (means = [1.0, 0.0],
                                                        vars = [1.0, 1.0],
                                                        weights = [0.5, 0.5])))
        @test_throws ArgumentError cluster(smld_dummy,
            MRFDensityClusterConfig(regime_gaussians = (means = [0.0, 1.0],
                                                        vars = [1.0, -1.0],
                                                        weights = [0.5, 0.5])))
        @test_throws ArgumentError cluster(smld_dummy,
            MRFDensityClusterConfig(regime_gaussians = (means = [0.0, 1.0],
                                                        vars = [1.0, 1.0],
                                                        weights = [0.5, 0.0])))
    end

    @testset "use_3d=true raises ArgumentError" begin
        smld = _make_2d_smld([(0.0, 0.0, 1), (1.0, 0.0, 1), (0.0, 1.0, 1)])
        @test_throws ArgumentError cluster(smld, MRFDensityClusterConfig(use_3d = true))
    end

    @testset "duplicate (x,y) coordinates raise ArgumentError" begin
        # Need ≥3 points per group to bypass the small-group short-circuit.
        pts = [(0.0, 0.0, 1), (1.0, 0.0, 1), (0.0, 1.0, 1), (1.0, 1.0, 1),
               (0.0, 0.0, 1)]  # duplicate of first point
        smld = _make_2d_smld(pts; n_datasets = 1)
        @test_throws ArgumentError cluster(smld, MRFDensityClusterConfig(per_dataset = false))
    end

    @testset "empty SMLD: zero clusters, empty regime vector" begin
        smld = _make_2d_smld(Tuple{Float64,Float64,Int}[]; n_datasets = 1)
        smld_out, info = cluster(smld, MRFDensityClusterConfig())
        @test info.n_locs_in == 0
        @test info.n_clustered == 0
        @test info.n_clusters == 0
        @test info.algorithm === :mrf_density
        @test smld_out.metadata["mrf_regime_per_emitter"] == Int[]
    end

    @testset "non-mutation: input emitter ids untouched" begin
        rng = Xoshiro(7)
        pts = Tuple{Float64,Float64,Int}[]
        for _ in 1:30
            push!(pts, (2.0 * rand(rng), 2.0 * rand(rng), 1))
        end
        # Add a tighter blob.
        append!(pts, _blob(rng, 1.0, 1.0, 0.05, 30))
        smld = _make_2d_smld(pts; n_datasets = 1)
        @test all(e -> e.id == 0, smld.emitters)
        cluster(smld, MRFDensityClusterConfig(per_dataset = false))
        @test all(e -> e.id == 0, smld.emitters)  # input unchanged
    end

    @testset "exports" begin
        @test :MRFDensityClusterConfig in names(SMLMClustering)
        @test :calibrate_regime_gaussians in names(SMLMClustering)
        @test :calibrate_regime_thresholds in names(SMLMClustering)
    end

    @testset "soft Gaussian unaries allow neighbor rescue" begin
        x = [0.49, 1.0, 1.0]
        thresholds = [0.5]
        gaussians = (means = [0.0, 1.0],
                     vars = [1.0, 1.0],
                     weights = [0.5, 0.5])
        adj = [[2, 3], [1, 3], [1, 2]]

        labels_hard = zeros(Int, 3)
        U_hard = SMLMClustering._unary_from_thresholds(x, thresholds, 2)
        SMLMClustering._icm_potts!(labels_hard, U_hard, adj, 1.0, 10)
        @test labels_hard == [1, 2, 2]

        labels_soft = zeros(Int, 3)
        U_soft = SMLMClustering._unary_from_gmm(x, gaussians)
        SMLMClustering._icm_potts!(labels_soft, U_soft, adj, 1.0, 10)
        @test labels_soft == [2, 2, 2]
    end

    @testset "Bayes boundary: weighted densities equal at the threshold" begin
        lwn(x, m, v, w) = log(w) - 0.5 * log(2π * v) - (x - m)^2 / (2v)
        bnd = SMLMClustering._gaussian_decision_boundary
        # Hand case: means 0 and 3, variances 1 and 0.25, equal weights. The root of
        # -3x² + 24x - 36 + log 4 = 0 in [0, 3] is 1.8877 (the old sign gave 2.119).
        t = bnd(0.0, 1.0, 0.5, 3.0, 0.25, 0.5)
        @test t ≈ (24 - sqrt(576 - 12 * (36 - log(4)))) / 6 atol = 1e-12
        @test lwn(t, 0.0, 1.0, 0.5) ≈ lwn(t, 3.0, 0.25, 0.5) atol = 1e-10
        # Swapped variances mirror it about 1.5.
        @test bnd(0.0, 0.25, 0.5, 3.0, 1.0, 0.5) ≈ 3.0 - t atol = 1e-12
        # Random two-component fits: wherever the weighted densities cross between
        # the means, the returned threshold is that crossing (posterior 0.5).
        rng = Xoshiro(11)
        errs = Float64[]; in_range = Bool[]
        for _ in 1:500
            m1 = 4 * rand(rng); m2 = m1 + 0.5 + 3 * rand(rng)
            v1 = 0.05 + 2 * rand(rng); v2 = 0.05 + 2 * rand(rng)
            w1 = 0.1 + 0.8 * rand(rng); w2 = 1 - w1
            f = x -> lwn(x, m1, v1, w1) - lwn(x, m2, v2, w2)
            f(m1) > 0 && f(m2) < 0 || continue             # one crossing lies in [m1, m2]
            x = bnd(m1, v1, w1, m2, v2, w2)
            push!(in_range, m1 <= x <= m2)
            push!(errs, abs(f(x)))
        end
        @test length(errs) > 300
        @test all(in_range)
        @test maximum(errs) < 1e-8
    end

    @testset "calibrate_regime_thresholds: argument validation" begin
        smld = _make_2d_smld([(0.0, 0.0, 1), (1.0, 0.0, 1), (0.0, 1.0, 1),
                              (1.0, 1.0, 1), (0.5, 0.5, 1)])
        @test_throws ArgumentError calibrate_regime_thresholds(smld; n_regimes = 1)
        @test_throws ArgumentError calibrate_regime_thresholds(smld; density_estimator = :bogus)
        @test_throws ArgumentError calibrate_regime_thresholds(smld; density_k = 0)
        # Too few emitters for n_regimes.
        smld_tiny = _make_2d_smld([(0.0, 0.0, 1)])
        @test_throws ArgumentError calibrate_regime_thresholds(smld_tiny; n_regimes = 2)
        @test_throws ArgumentError calibrate_regime_gaussians(smld; n_regimes = 1)
        @test_throws ArgumentError calibrate_regime_gaussians(smld; density_estimator = :bogus)
        @test_throws ArgumentError calibrate_regime_gaussians(smld; density_k = 0)
        @test_throws ArgumentError calibrate_regime_gaussians(smld_tiny; n_regimes = 2)
    end

end
