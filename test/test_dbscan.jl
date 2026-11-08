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

# Independent brute-force classical DBSCAN (Ester et al. 1996) for a 2 × n matrix:
# distance-matrix neighborhoods (d ≤ r, self included), core iff |N| ≥ minpts, core
# components by BFS in index order, a border point to the lowest core-cluster id
# within r, clusters left below minpts dropped and the rest renumbered.
# Returns (labels, core).
function _ref_dbscan_classical(X::AbstractMatrix, r::Real, minpts::Int)
    n = size(X, 2)
    N = [[j for j in 1:n if sqrt((X[1, i] - X[1, j])^2 + (X[2, i] - X[2, j])^2) <= r]
         for i in 1:n]
    core = [length(N[i]) >= minpts for i in 1:n]
    lab = zeros(Int, n); k = 0
    for i in 1:n
        (core[i] && lab[i] == 0) || continue
        k += 1; lab[i] = k; queue = [i]
        while !isempty(queue)
            p = popfirst!(queue)
            for q in N[p]
                core[q] && lab[q] == 0 && (lab[q] = k; push!(queue, q))
            end
        end
    end
    for i in 1:n
        core[i] && continue
        cids = [lab[j] for j in N[i] if core[j]]
        lab[i] = isempty(cids) ? 0 : minimum(cids)
    end
    sizes = [count(==(c), lab) for c in 1:k]
    newid = cumsum(sizes .>= minpts)
    return [l == 0 || sizes[l] < minpts ? 0 : newid[l] for l in lab], core
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

    if SMLM_TEST_FULL
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
    end  # SMLM_TEST_FULL

    if SMLM_TEST_FULL
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
    end  # SMLM_TEST_FULL

    if SMLM_TEST_FULL
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
    end  # SMLM_TEST_FULL

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

    if SMLM_TEST_FULL
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
    end  # SMLM_TEST_FULL

    if SMLM_TEST_FULL
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
    end  # SMLM_TEST_FULL

end
