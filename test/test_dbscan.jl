using SMLMClustering
using SMLMData
using Test
using Random
include(joinpath(@__DIR__, "long", "utils", "helpers.jl"))

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

    @testset "a border point's own neighbors are not added (classical DBSCAN)" begin
        # eps = 100 nm, min_points = 4. P1..P4 form a core blob; Q is 90 nm from P3
        # and > 100 nm from P1, P2, P4, so N(Q) = {Q, P3, R} (3 < 4): Q is a border
        # point. R is 90 nm from Q and 180 nm from P3: within eps of no core point,
        # so classical DBSCAN labels it noise. Clustering.jl 0.15's dbscan put R in
        # the cluster (sizes [6]).
        P1 = (0.980, 1.000, 1); P2 = (0.985, 1.020, 1); P3 = (1.000, 1.000, 1)
        P4 = (0.985, 0.980, 1); Q = (1.090, 1.000, 1); R = (1.180, 1.000, 1)
        cfg = DBSCANConfig(eps_nm = 100.0, min_points = 4, per_dataset = false)
        out, info = cluster(_make_2d_smld([P1, P2, P3, P4, Q, R]), cfg)
        @test [e.id for e in out.emitters] == [1, 1, 1, 1, 1, 0]
        @test info.cluster_sizes == [5]
        @test info.n_noise == 1
        # Same answer whatever the emitter order.
        out2, info2 = cluster(_make_2d_smld([R, Q, P4, P3, P2, P1]), cfg)
        @test [e.id for e in out2.emitters] == [0, 1, 1, 1, 1, 1]
        @test info2.cluster_sizes == [5]
    end

    @testset "matches a brute-force classical DBSCAN on random data" begin
        ref_dbscan = _ref_dbscan_classical
        rng = Xoshiro(20261107)
        pts = Tuple{Float64,Float64,Int}[]
        for c in 1:6
            append!(pts, _blob(rng, 0.4 + 0.6 * rand(rng), 0.4 + 0.6 * rand(rng), 0.04, 40))
        end
        for _ in 1:150
            push!(pts, (1.4 * rand(rng), 1.4 * rand(rng), 1))
        end
        smld = _make_2d_smld(pts)
        X = reduce(hcat, [[p[1], p[2]] for p in pts])          # 2 × n, μm
        n = size(X, 2)
        D = [hypot(X[1, i] - X[1, j], X[2, i] - X[2, j]) for i in 1:n, j in 1:n]
        n_halo = 0
        for (eps, mp) in ((30.0, 2), (40.0, 3), (40.0, 5), (50.0, 5), (60.0, 8))
            out, info = cluster(smld, DBSCANConfig(eps_nm = eps, min_points = mp, per_dataset = false))
            ids = [e.id for e in out.emitters]
            ref, core = ref_dbscan(X, eps / 1000, mp)
            @test ids == ref
            @test info.n_clustered == count(!=(0), ref)
            @test info.cluster_sizes == [count(==(c), ref) for c in 1:maximum(ref; init = 0)]
            # Every clustered non-core point lies within eps of a core point of its cluster.
            r = eps / 1000
            @test all(i -> ids[i] == 0 || core[i] ||
                           any(j -> core[j] && ids[j] == ids[i] && D[i, j] <= r, 1:n), 1:n)
            # Noise points within eps of a border point that is itself within eps of a
            # core point: the points a border-of-border expansion would wrongly add.
            border = [!core[j] && any(q -> core[q] && D[j, q] <= r, 1:n) for j in 1:n]
            n_halo += count(i -> ref[i] == 0 && !core[i] &&
                                 any(j -> border[j] && D[i, j] <= r, 1:n), 1:n)
        end
        @test n_halo > 0        # the data set exercises the old bug's path
    end

    @testset "agrees with PrecisionDBSCANConfig at constant precision" begin
        # σ_x = σ_y = 10 nm for every emitter, so PrecisionDBSCAN's threshold
        # nsigma·(σ_i + σ_j) is the constant 2.5·20 nm = 50 nm (strict <); DBSCAN
        # uses ≤ 50 nm. With continuous random coordinates no pair sits at 50 nm.
        rng = Xoshiro(7)
        pts = Tuple{Float64,Float64,Int}[]
        for _ in 1:5
            append!(pts, _blob(rng, 2.0 * rand(rng), 2.0 * rand(rng), 0.05, 50))
        end
        for _ in 1:200
            push!(pts, (2.0 * rand(rng), 2.0 * rand(rng), 1))
        end
        smld = _make_2d_smld(pts)
        for mp in (3, 5, 8)
            a, ia = cluster(smld, DBSCANConfig(eps_nm = 50.0, min_points = mp, per_dataset = false))
            b, ib = cluster(smld, PrecisionDBSCANConfig(nsigma = 2.5, min_points = mp, per_dataset = false))
            @test [e.id for e in a.emitters] == [e.id for e in b.emitters]
            @test ia.cluster_sizes == ib.cluster_sizes
        end
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

end
