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

    @testset "Bayes boundary outside the two means" begin
        lwn(x, m, v, w) = log(w) - 0.5 * log(2π * v) - (x - m)^2 / (2v)
        bnd = SMLMClustering._gaussian_decision_boundary
        # Means 0 and 3, variances 1 and 0.25, the narrow upper component 400 times
        # lighter: it wins only on [3.32, 4.68], beyond its own mean, so the boundary is
        # 3.32, where it overtakes. The released code returned the midpoint 1.5.
        w1, w2 = 400 / 401, 1 / 401
        C = -36 + log(4) - 2 * log(w1 / w2)
        t = bnd(0.0, 1.0, w1, 3.0, 0.25, w2)
        @test t ≈ (24 - sqrt(576 + 12 * C)) / 6 atol = 1e-12
        @test t > 3.0
        @test lwn(t, 0.0, 1.0, w1) ≈ lwn(t, 3.0, 0.25, w2) atol = 1e-10
        # Mirror: a narrow, light lower component; the boundary falls below the lower mean.
        @test bnd(0.0, 0.25, w2, 3.0, 1.0, w1) ≈ 3.0 - t atol = 1e-12
        # Random fits, those with the crossing outside the means included: wherever the
        # weighted densities cross, the threshold is the crossing where the upper component
        # overtakes the lower one.
        rng = Xoshiro(12)
        errs = Float64[]; rising = Bool[]; n_outside = 0
        for _ in 1:2000
            m1 = 4 * rand(rng); m2 = m1 + 0.5 + 3 * rand(rng)
            v1 = 0.05 + 2 * rand(rng); v2 = 0.05 + 2 * rand(rng)
            a = 0.02 + 0.96 * rand(rng); b = 1 - a
            f = x -> lwn(x, m1, v1, a) - lwn(x, m2, v2, b)       # > 0: the lower one wins
            xs = range(m1 - 30, m2 + 30; length = 20001)
            (any(x -> f(x) > 0, xs) && any(x -> f(x) < 0, xs)) || continue
            n_outside += (f(m1) > 0) == (f(m2) > 0)
            x = bnd(m1, v1, a, m2, v2, b)
            push!(errs, abs(f(x)))
            push!(rising, f(x - 1e-6) > f(x + 1e-6))
        end
        @test n_outside > 100
        @test maximum(errs) < 1e-8
        @test all(rising)
    end

    @testset "regime thresholds: a component that never wins gets an empty regime" begin
        lwn(x, m, v, w) = log(w) - 0.5 * log(2π * v) - (x - m)^2 / (2v)
        thr = SMLMClustering._regime_thresholds
        # Two components: the pairwise boundary.
        @test thr((means = [0.0, 3.0], vars = [1.0, 0.25], weights = [0.5, 0.5])) ==
              [SMLMClustering._gaussian_decision_boundary(0.0, 1.0, 0.5, 3.0, 0.25, 0.5)]
        # Three equal-variance components, the middle one 50 times lighter: the pairwise
        # boundaries 1|2 = 2.95 and 2|3 = 1.05 are out of order (the released code returned
        # them, and cluster rejects unsorted thresholds), so component 2 never wins and both
        # thresholds are the 1|3 boundary, 2.0.
        fit3 = (means = [0.0, 2.0, 4.0], vars = [1.0, 1.0, 1.0], weights = [0.495, 0.01, 0.495])
        @test thr(fit3) ≈ [2.0, 2.0] atol = 1e-12
        # Equal variances: the thresholds reproduce the Bayes classifier argmax_k w_k N_k
        # exactly (bins as in _unary_from_thresholds: x >= t_k is above regime k).
        rng = Xoshiro(13)
        agree = Bool[]; sorted = Bool[]; n_empty = 0
        for _ in 1:300
            K = rand(rng, 3:5)
            m = sort(6 * rand(rng, K)); v = fill(0.1 + rand(rng), K)
            w = 0.02 .+ rand(rng, K); w ./= sum(w)
            t = thr((means = m, vars = v, weights = w))
            push!(sorted, issorted(t))
            n_empty += any(k -> t[k] == t[k+1], 1:K-2)
            for x in range(m[1] - 3, m[end] + 3; length = 401)
                bayes = argmax([lwn(x, m[k], v[k], w[k]) for k in 1:K])
                push!(agree, 1 + count(<=(x), t) == bayes)
            end
        end
        @test n_empty > 20
        @test all(agree)
        # Unequal variances: always sorted, so cluster accepts them.
        for _ in 1:300
            K = rand(rng, 3:5)
            m = sort(6 * rand(rng, K)); v = 0.05 .+ 2 * rand(rng, K)
            w = 0.02 .+ rand(rng, K); w ./= sum(w)
            push!(sorted, issorted(thr((means = m, vars = v, weights = w))))
        end
        @test all(sorted)
    end

    @testset "regime thresholds: no crossing, one component wins everywhere" begin
        lwn(x, m, v, w) = log(w) - 0.5 * log(2π * v) - (x - m)^2 / (2v)
        bnd = SMLMClustering._gaussian_decision_boundary
        thr = SMLMClustering._regime_thresholds
        # The threshold is where the upper regime starts: -Inf when the upper component
        # wins everywhere, +Inf when it never does. Same mean and variance: the heavier wins
        # everywhere; identical components: no preference, any split (the shared mean).
        @test bnd(1.0, 1.0, 0.3, 1.0, 1.0, 0.7) == -Inf
        @test bnd(1.0, 1.0, 0.7, 1.0, 1.0, 0.3) == Inf
        @test bnd(1.0, 1.0, 0.5, 1.0, 1.0, 0.5) == 1.0
        # Unequal variances, no crossing: a narrow, light component under a broad one.
        f = x -> lwn(x, 0.0, 1.0, 0.999) - lwn(x, 1.0, 0.01, 0.001)
        @test all(x -> f(x) > 0, range(-50, 50; length = 10001))
        @test bnd(0.0, 1.0, 0.999, 1.0, 0.01, 0.001) == Inf
        @test bnd(-1.0, 0.01, 0.001, 0.0, 1.0, 0.999) == -Inf
        # Nearly equal variances (A ≈ 1e-13, not zero): the densities cross far out, and
        # the threshold is that crossing, not the midpoint. The inputs are ill-conditioned
        # (A loses digits in Float64), hence the loose tolerance against a BigFloat root.
        p = (0.0, 1.0, 0.9, 1e-14, 1 + 1e-13, 0.1)
        A = 1 / big(p[2]) - 1 / big(p[5])
        B = -2 * (big(p[1]) / p[2] - big(p[4]) / p[5])
        C = big(p[1])^2 / p[2] - big(p[4])^2 / p[5] + log(big(p[2]) / p[5]) - 2 * log(big(p[3]) / p[6])
        @test bnd(p...) ≈ Float64((-B + sqrt(B^2 - 4A * C)) / (2A)) rtol = 1e-2
        @test bnd(p...) > 1e6
        # A narrow, light middle component under two broad ones never wins: its regime is
        # empty and both thresholds are the 1|3 boundary, 1.0 (not the midpoints 0.5, 1.5).
        @test thr((means = [0.0, 1.0, 2.0], vars = [1.0, 0.01, 1.0],
                   weights = [0.4995, 0.001, 0.4995])) ≈ [1.0, 1.0] atol = 1e-12
        # Same mean and variance as its lower neighbor, lighter: component 3 never wins.
        t = thr((means = [0.0, 1.0, 1.0, 3.0], vars = ones(4), weights = [0.1, 0.6, 0.2, 0.1]))
        @test t[1] ≈ bnd(0.0, 1.0, 0.1, 1.0, 1.0, 0.6) atol = 1e-12
        @test t[2] == t[3]
        @test t[2] ≈ bnd(1.0, 1.0, 0.6, 3.0, 1.0, 0.1) atol = 1e-12
        # An end component that never wins: its regime is empty through an infinite threshold.
        @test thr((means = [0.0, 1.0], vars = [1.0, 0.01], weights = [0.999, 0.001])) == [Inf]
        t = thr((means = [-1.0, 0.0, 3.0], vars = [0.01, 1.0, 1.0], weights = [0.001, 0.5, 0.499]))
        @test t[1] == -Inf
        @test t[2] ≈ bnd(0.0, 1.0, 0.5, 3.0, 1.0, 0.499) atol = 1e-12
        # Dominated by a component that is not its neighbor: the broad first component is
        # above the narrow, light third one everywhere (pair boundary +Inf), so component 3
        # gets no regime, though it overtakes component 2 at 2.47. Mirrored: component 1.
        @test bnd(0.0, 1.0, 0.5, 2.0, 0.25, 0.01) == Inf
        t = thr((means = [0.0, 1.0, 2.0], vars = [1.0, 0.25, 0.25], weights = [0.5, 0.49, 0.01]))
        @test t[1] ≈ bnd(0.0, 1.0, 0.5, 1.0, 0.25, 0.49) atol = 1e-12
        @test t[2] == Inf
        t = thr((means = [0.0, 1.0, 2.0], vars = [0.25, 0.25, 1.0], weights = [0.01, 0.49, 0.5]))
        @test t[1] == -Inf
        @test t[2] ≈ bnd(1.0, 0.25, 0.49, 2.0, 1.0, 0.5) atol = 1e-12
        # Unequal variances: every component that some other component is above everywhere
        # gets an empty regime, and the thresholds stay sorted. Oracle: lwn_j - lwn_k is the
        # quadratic a x² + b x + c below; j is above k everywhere iff a > 0 and b² < 4ac.
        function above_everywhere(j, k, m, v, w)
            a = 1 / (2v[k]) - 1 / (2v[j])
            b = m[j] / v[j] - m[k] / v[k]
            c = log(w[j] / w[k]) - log(v[j] / v[k]) / 2 - m[j]^2 / (2v[j]) + m[k]^2 / (2v[k])
            return a > 0 && b^2 < 4a * c
        end
        rng = Xoshiro(15)
        ok = Bool[]; n_dom = 0
        for _ in 1:500
            K = rand(rng, 2:5)
            m = sort(4 * rand(rng, K)); v = 0.05 .+ 2 * rand(rng, K)
            w = 0.005 .+ rand(rng, K); w ./= sum(w)
            t = thr((means = m, vars = v, weights = w))
            te = [-Inf; t; Inf]                    # regime k is [te[k], te[k+1])
            for k in 1:K
                dominated = any(j -> j != k && above_everywhere(j, k, m, v, w), 1:K)
                n_dom += dominated
                dominated && push!(ok, te[k] == te[k + 1])
            end
            push!(ok, issorted(t))
        end
        @test n_dom > 50
        @test all(ok)
        # Equal variances with tied means: the thresholds still reproduce the Bayes
        # classifier argmax_k w_k N_k (bins as in _unary_from_thresholds).
        rng = Xoshiro(14)
        agree = Bool[]; n_tied = 0
        for _ in 1:300
            K = rand(rng, 2:5)
            m = sort(Float64.(rand(rng, 0:3, K))); v = fill(0.1 + rand(rng), K)
            w = 0.02 .+ rand(rng, K); w ./= sum(w)
            n_tied += !allunique(m)
            t = thr((means = m, vars = v, weights = w))
            for x in range(m[1] - 3, m[end] + 3; length = 401)
                bayes = argmax([lwn(x, m[k], v[k], w[k]) for k in 1:K])
                push!(agree, 1 + count(<=(x), t) == bayes)
            end
        end
        @test n_tied > 100
        @test all(agree)
    end

    @testset "regime thresholds: every non-empty regime holds points its component wins" begin
        lwn(x, m, v, w) = log(w) - 0.5 * log(2π * v) - (x - m)^2 / (2v)
        thr = SMLMClustering._regime_thresholds
        # Brute force: each non-empty regime [t_{k-1}, t_k) holds a point where component k
        # has the highest weighted density of all the components (ties count). The points
        # are a fine grid plus one probe inside every regime, so a regime that lies off the
        # grid is still checked.
        function wins_in_own_regime(m, v, w, t)
            K = length(m); te = [-Inf; t; Inf]
            top(x) = maximum(j -> lwn(x, m[j], v[j], w[j]), 1:K)
            grid = collect(range(minimum(m .- 12 .* sqrt.(v)), maximum(m .+ 12 .* sqrt.(v)); length = 20001))
            gtop = map(top, grid)
            return all(1:K) do k
                a, b = te[k], te[k + 1]
                a < b || return true                                   # empty regime
                probe = isfinite(a) && isfinite(b) ? (a + b) / 2 :
                        isfinite(a) ? a + 1 + abs(a) : isfinite(b) ? b - 1 - abs(b) : 0.0
                any(i -> a <= grid[i] < b && lwn(grid[i], m[k], v[k], w[k]) >= gtop[i], eachindex(grid)) ||
                    lwn(probe, m[k], v[k], w[k]) >= top(probe)
            end
        end
        # The oracle itself: a regime entirely off the grid is still checked (component 1 is
        # below component 2 everywhere, so a cut at -100 must fail), and a tie between
        # identical components counts as winning.
        mT, vT, wT = [0.0, 1.0], [0.01, 1.0], [0.001, 0.999]
        @test !wins_in_own_regime(mT, vT, wT, [-100.0])
        @test thr((means = mT, vars = vT, weights = wT)) == [-Inf]
        fI = (means = [1.0, 1.0], vars = [1.0, 1.0], weights = [0.5, 0.5])
        @test wins_in_own_regime(fI.means, fI.vars, fI.weights, thr(fI))
        # Touching the top at a single point is not winning: the narrow component equals
        # the broad one only at 0 and is below it everywhere else, so it gets no regime.
        @test thr((means = [0.0, 0.0], vars = [1.0, 4.0], weights = [1.0, 2.0])) == [-Inf]
        # 4 regimes, unequal variances: the neighbor-only walk gave component 3 the range
        # [-0.83, -0.21), where component 1 is highest everywhere; now its regime is empty.
        w4 = [0.342, 0.254, 0.099, 0.304]; w4 ./= sum(w4)
        f4 = (means = [-2.506, -2.313, 0.291, 1.650], vars = [2.447, 0.496, 4.875, 0.858], weights = w4)
        t4 = thr(f4)
        @test t4[2] == t4[3]
        @test wins_in_own_regime(f4.means, f4.vars, f4.weights, t4)
        # 3 regimes: the end component 1 kept (-Inf, -0.83) though it is highest nowhere,
        # while the broad component 2, highest on most of that range, had an empty regime.
        f3 = (means = [0.433, 0.515, 0.729], vars = [0.826, 1.802, 0.663], weights = [0.234, 0.267, 0.499])
        t3 = thr(f3)
        @test t3[1] == -Inf
        @test -0.71 < t3[2] < -0.69
        @test wins_in_own_regime(f3.means, f3.vars, f3.weights, t3)
        # Random unequal-variance mixtures at 3 to 5 regimes.
        rng = Xoshiro(16)
        ok = Bool[]
        for K in 3:5, _ in 1:400
            m = sort(2 .* randn(rng, K)); v = exp.(randn(rng, K)); w = rand(rng, K); w ./= sum(w)
            t = thr((means = m, vars = v, weights = w))
            push!(ok, issorted(t) && wins_in_own_regime(m, v, w, t))
        end
        @test all(ok)
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
