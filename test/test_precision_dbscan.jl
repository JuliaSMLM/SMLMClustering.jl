using SMLMClustering
# The reuse primitive is public but not exported — bring the names into scope explicitly.
using SMLMClustering: build_precision_neighbor_graph, precision_dbscan_labels,
                      precision_dbscan_labels!, PrecisionNeighborGraph
using SMLMData
using Clustering        # to cross-check min_points against classical DBSCAN minPts
using Test
using Random

# Reuse the 2D SMLD builder shape from test_dbscan.jl (precision σ_x = σ_y = 0.01 μm).
function _make_prec_smld(points::Vector{Tuple{Float64,Float64,Int}};
                         σ::Float64 = 0.01,
                         n_datasets::Int = maximum(p[3] for p in points))
    cam = IdealCamera(1:64, 1:64, 0.1)
    emitters = [Emitter2DFit{Float64}(
        x, y, 1000.0, 10.0, σ, σ, 50.0, 2.0;
        frame = 1, dataset = ds,
    ) for (x, y, ds) in points]
    BasicSMLD(emitters, cam, 1, n_datasets, Dict{String,Any}())
end

_prec_blob(rng, cx, cy, σ, n; dataset = 1) =
    [(cx + σ * randn(rng), cy + σ * randn(rng), dataset) for _ in 1:n]

@testset "Precision-weighted DBSCAN" begin

    # ---- primitive: neighbor graph -----------------------------------------
    @testset "build_precision_neighbor_graph — CSR + geometry" begin
        coords = [0.0 1.0 2.0; 0.0 0.0 0.0]      # 3 collinear points, spacing 1
        g = build_precision_neighbor_graph(coords, 2.5)
        @test g isa PrecisionNeighborGraph
        @test g.n == 3
        @test g.dims == 2
        @test g.max_radius == 2.5
        @test g.offsets == [1, 3, 5, 7]           # full adjacency within 2.5
        @test g.neighbors == [2, 3, 1, 3, 1, 2]   # index-sorted per point
        @test g.dists ≈ [1.0, 2.0, 1.0, 1.0, 2.0, 1.0]

        # tighter radius drops the far pair (1↔3, d=2)
        g2 = build_precision_neighbor_graph(coords, 1.5)
        @test g2.offsets == [1, 2, 4, 5]
        @test g2.neighbors == [2, 1, 3, 2]

        # empty input is well-formed
        ge = build_precision_neighbor_graph(zeros(2, 0), 1.0)
        @test ge.n == 0 && isempty(ge.neighbors) && ge.offsets == [1]

        @test_throws ArgumentError build_precision_neighbor_graph(zeros(2, 3), 0.0)
        @test_throws ArgumentError build_precision_neighbor_graph(zeros(4, 3), 1.0)  # bad dims
    end

    # ---- primitive: label pass ---------------------------------------------
    @testset "precision_dbscan_labels — thresholds + branches" begin
        coords = [0.0 1.0 2.0; 0.0 0.0 0.0]
        g = build_precision_neighbor_graph(coords, 2.5)
        σ = [1.0, 1.0, 1.0]

        # nsigma=0.6 → active iff d < 0.6*(1+1)=1.2 → the two unit edges → one CC
        @test precision_dbscan_labels(g, σ, 0.6; min_points = 0) == [1, 1, 1]
        # nsigma=0.4 → threshold 0.8 → no active edge → three singletons
        @test precision_dbscan_labels(g, σ, 0.4; min_points = 0) == [1, 2, 3]
        # core-point, self-inclusive minPts: min_points=2 → core iff active degree ≥ 1,
        # so all three points are core → one cluster
        @test precision_dbscan_labels(g, σ, 0.6; min_points = 2) == [1, 1, 1]

        # add an isolated 4th point → noise under the core-point branch
        coords4 = [0.0 1.0 2.0 10.0; 0.0 0.0 0.0 0.0]
        g4 = build_precision_neighbor_graph(coords4, 2.5)
        @test precision_dbscan_labels(g4, fill(1.0, 4), 0.6; min_points = 2) == [1, 1, 1, 0]
        # min_points=0 labels every point (singleton gets its own id)
        @test precision_dbscan_labels(g4, fill(1.0, 4), 0.6; min_points = 0) == [1, 1, 1, 2]

        # in-place form matches the allocating form
        buf = Vector{Int}(undef, 3)
        @test precision_dbscan_labels!(buf, g, σ, 0.6; min_points = 0) === buf
        @test buf == [1, 1, 1]

        # superset guard: nsigma*2*max(σ) must not exceed max_radius
        @test_throws ArgumentError precision_dbscan_labels(g, σ, 2.0; min_points = 0)  # need 4.0 > 2.5
        @test_throws ArgumentError precision_dbscan_labels(g, [1.0, 1.0], 0.6)      # σ length ≠ n
        @test_throws ArgumentError precision_dbscan_labels(g, σ, 0.0)               # nsigma ≤ 0
    end

    # ---- the reuse invariant that BaGoL depends on -------------------------
    @testset "cache reuse: coarse-built graph == tight-built graph" begin
        rng = Xoshiro(7)
        pts = [(0.10 + 0.02 * randn(rng), 0.10 + 0.02 * randn(rng)) for _ in 1:40]
        append!(pts, [(0.60 + 0.02 * randn(rng), 0.60 + 0.02 * randn(rng)) for _ in 1:40])
        coords = Matrix(reshape(reduce(vcat, ([p[1], p[2]] for p in pts)), 2, :))  # 2×N
        n = size(coords, 2)
        σ = fill(0.03, n)
        for ns in (2.0, 3.0, 4.0), mp in (0, 3)
            tight = build_precision_neighbor_graph(coords, ns * 2 * maximum(σ))
            coarse = build_precision_neighbor_graph(coords, 5.0 * 2 * maximum(σ))  # superset
            lt = precision_dbscan_labels(tight, σ, ns; min_points = mp)
            lc = precision_dbscan_labels(coarse, σ, ns; min_points = mp, check_superset = false)
            @test lt == lc                       # identical partition AND canonical ids
        end
    end

    @testset "strict boundary + heterogeneous-σ reuse + zero σ" begin
        # d == threshold is NOT active (strict `<`): two points 1.0 apart, σ=0.5 each,
        # nsigma=1.0 → threshold exactly 1.0 → not linked.
        gb = build_precision_neighbor_graph([0.0 1.0; 0.0 0.0], 1.0)
        @test precision_dbscan_labels(gb, [0.5, 0.5], 1.0; min_points = 0) == [1, 2]
        gb2 = build_precision_neighbor_graph([0.0 1.0; 0.0 0.0], 1.5)
        @test precision_dbscan_labels(gb2, [0.5, 0.5], 1.01; min_points = 0) == [1, 1]  # d=1.0 < 1.01

        # reuse invariant holds with per-point heterogeneous σ_eff (the realistic case)
        rng = Xoshiro(99)
        coords = Matrix(reshape(reduce(vcat,
            ([0.1 + 0.03 * randn(rng), 0.1 + 0.03 * randn(rng)] for _ in 1:60)), 2, :))
        n = size(coords, 2)
        σ = [0.02 + 0.03 * (i / n) for i in 1:n]
        for ns in (2.0, 3.5), mp in (0, 3)
            tight = build_precision_neighbor_graph(coords, ns * 2 * maximum(σ))
            coarse = build_precision_neighbor_graph(coords, 4.0 * 2 * maximum(σ))
            @test precision_dbscan_labels(tight, σ, ns; min_points = mp) ==
                  precision_dbscan_labels(coarse, σ, ns; min_points = mp, check_superset = false)
        end

        # all-zero precision → config cannot form a neighborhood → ArgumentError
        smld0 = _make_prec_smld([(0.0, 0.0, 1), (0.01, 0.0, 1)]; σ = 0.0)
        @test_throws ArgumentError cluster(smld0, PrecisionDBSCANConfig(nsigma = 5.0, per_dataset = false))
    end

    @testset "core-point branch: multi-cluster ids + border tie-break" begin
        # two separated pairs → two components; min_points=1 → all core; canonical ids by index
        gc = build_precision_neighbor_graph([0.0 0.01 1.0 1.01; 0.0 0.0 0.0 0.0], 0.2)
        @test precision_dbscan_labels(gc, fill(0.1, 4), 1.0; min_points = 1) == [1, 1, 2, 2]
        @test precision_dbscan_labels(gc, fill(0.1, 4), 1.0; min_points = 0) == [1, 1, 2, 2]

        # Border adjacent to two distinct core clusters joins the LOWER-id one.
        # A = idx 1-5 (blob at x≈0 + antenna idx5), B = idx 6-10 (blob at x≈1 + antenna
        # idx10). P = idx11 midway with large σ (0.40) reaches one antenna in each cluster,
        # active degree 2. With self-inclusive min_points=4 (core iff active degree ≥ 3),
        # P (degree 2) is a *border*, not a core bridge, so A and B stay disjoint. A is
        # labeled 1 (first core by index), B is 2 → P joins 1.
        xs = [0.0, 0.02, 0.0, 0.02, 0.07, 1.0, 1.02, 1.0, 1.02, 0.93, 0.5]
        ys = [0.0, 0.0, 0.02, 0.02, 0.01, 0.0, 0.0, 0.02, 0.02, 0.01, 0.01]
        coords = permutedims(hcat(xs, ys))            # 2×11, columns = points
        σ = [fill(0.05, 10); 0.40]
        gp = build_precision_neighbor_graph(coords, 0.80)
        @test precision_dbscan_labels(gp, σ, 1.0; min_points = 4) ==
              [1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 1]
    end

    @testset "min_points == classical DBSCAN minPts (self-inclusive)" begin
        # 5 collinear points (radius 1.5 → interior active degree 2, endpoints 1) + 1
        # isolated point. The self-inclusive fix makes precision's min_points identical
        # to Clustering.dbscan's min_neighbors across the core/border/noise regimes:
        # min_points=3 clusters the whole chain (it was ALL-NOISE under the old
        # self-exclusive count) and marks the isolated point noise.
        X = [0.0 1.0 2.0 3.0 4.0 20.0; 0.0 0.0 0.0 0.0 0.0 0.0]
        g = build_precision_neighbor_graph(X, 30.0)
        σ = fill(1.0, 6)                      # threshold = nsigma·2; nsigma=0.75 → 1.5
        for m in (2, 3, 4)
            ref = Clustering.dbscan(X, 1.5; min_neighbors = m, min_cluster_size = 1).assignments
            @test precision_dbscan_labels(g, σ, 0.75; min_points = m) == ref
        end
        @test precision_dbscan_labels(g, σ, 0.75; min_points = 3) == [1, 1, 1, 1, 1, 0]
    end

    # ---- SMLD-facing config ------------------------------------------------
    @testset "PrecisionDBSCANConfig construction + defaults" begin
        cfg = PrecisionDBSCANConfig(nsigma = 5.0)
        @test cfg isa AbstractClusterConfig
        @test cfg.nsigma == 5.0
        @test cfg.min_points == 5
        @test cfg.use_3d === false
        @test cfg.per_dataset === true
        @test cfg.remove_unclustered === false
    end

    @testset "config validation" begin
        smld = _make_prec_smld([(0.0, 0.0, 1), (1.0, 1.0, 1)])
        @test_throws ArgumentError cluster(smld, PrecisionDBSCANConfig(nsigma = 0.0))
        @test_throws ArgumentError cluster(smld, PrecisionDBSCANConfig(nsigma = 5.0, min_points = 0))
    end

    if SMLM_TEST_FULL
    @testset "two blobs recovered + non-mutating + emitter.id" begin
        rng = Xoshiro(20260702)
        pts = Tuple{Float64,Float64,Int}[]
        append!(pts, _prec_blob(rng, 1.0, 1.0, 0.005, 40))
        append!(pts, _prec_blob(rng, 5.0, 5.0, 0.005, 40))
        for k in 1:6
            push!(pts, (10.0 + k, 10.0 + k, 1))     # far, sparse noise
        end
        smld = _make_prec_smld(pts; σ = 0.01, n_datasets = 1)
        n_in = length(smld.emitters)

        cfg = PrecisionDBSCANConfig(nsigma = 5.0, min_points = 5, per_dataset = false)
        smld_out, info = cluster(smld, cfg)
        @test info isa ClusterInfo
        @test info.algorithm === :precision_dbscan
        @test info.n_locs_in == n_in
        @test info.n_clusters == 2
        @test sum(info.cluster_sizes) == info.n_clustered
        @test info.n_clustered + info.n_noise == info.n_locs_in
        @test info.n_clustered >= 78              # both 40-pt blobs cluster
        # non-mutating: input ids untouched, output carries labels
        @test all(e -> e.id == 0, smld.emitters)
        @test all(e -> e.id in 0:info.n_clusters, smld_out.emitters)
        for k in 1:info.n_clusters
            @test count(e -> e.id == k, smld_out.emitters) == info.cluster_sizes[k]
        end

        # remove_unclustered drops the noise
        _, info_rm = cluster(smld, PrecisionDBSCANConfig(nsigma = 5.0, min_points = 5,
                                                         per_dataset = false, remove_unclustered = true))
        smld_rm, _ = cluster(smld, PrecisionDBSCANConfig(nsigma = 5.0, min_points = 5,
                                                         per_dataset = false, remove_unclustered = true))
        @test all(e -> e.id != 0, smld_rm.emitters)
        @test length(smld_rm.emitters) == info_rm.n_clustered
    end
    end  # SMLM_TEST_FULL

    if SMLM_TEST_FULL
    @testset "per_dataset label namespace is local" begin
        rng = Xoshiro(11)
        pts = Tuple{Float64,Float64,Int}[]
        append!(pts, _prec_blob(rng, 1.0, 1.0, 0.005, 20; dataset = 1))
        append!(pts, _prec_blob(rng, 3.0, 3.0, 0.005, 20; dataset = 1))
        append!(pts, _prec_blob(rng, 1.0, 1.0, 0.005, 20; dataset = 2))
        append!(pts, _prec_blob(rng, 3.0, 3.0, 0.005, 20; dataset = 2))
        smld = _make_prec_smld(pts; σ = 0.01, n_datasets = 2)

        _, info = cluster(smld, PrecisionDBSCANConfig(nsigma = 5.0, min_points = 5, per_dataset = true))
        @test info.n_clusters == 4
        @test info.n_clustered == 80
        smld_out, _ = cluster(smld, PrecisionDBSCANConfig(nsigma = 5.0, min_points = 5, per_dataset = true))
        @test sort!(unique(e.id for e in smld_out.emitters if e.dataset == 1)) == [1, 2]
        @test sort!(unique(e.id for e in smld_out.emitters if e.dataset == 2)) == [1, 2]

        _, info_flat = cluster(smld, PrecisionDBSCANConfig(nsigma = 5.0, min_points = 5, per_dataset = false))
        @test info_flat.n_clusters == 2
        @test info_flat.n_clustered == 80
    end
    end  # SMLM_TEST_FULL

    if SMLM_TEST_FULL
    @testset "3D precision clustering + empty SMLD" begin
        cam = IdealCamera(1:64, 1:64, 0.1)
        rng = Xoshiro(3)
        emitters = SMLMData.Emitter3DFit{Float64}[]
        for zc in (0.0, 1.0), _ in 1:25
            push!(emitters, Emitter3DFit{Float64}(
                2.0 + 0.005 * randn(rng), 2.0 + 0.005 * randn(rng), zc + 0.005 * randn(rng),
                1000.0, 10.0, 0.01, 0.01, 0.01, 50.0, 2.0; frame = 1, dataset = 1))
        end
        smld3 = BasicSMLD(emitters, cam, 1, 1, Dict{String,Any}())
        _, info = cluster(smld3, PrecisionDBSCANConfig(nsigma = 5.0, min_points = 5,
                                                       use_3d = true, per_dataset = false))
        @test info.n_clusters == 2
        @test info.n_clustered == 50

        smld_empty = _make_prec_smld(Tuple{Float64,Float64,Int}[]; n_datasets = 1)
        smld_out, info_e = cluster(smld_empty, PrecisionDBSCANConfig(nsigma = 5.0, per_dataset = false))
        @test info_e.n_locs_in == 0 && info_e.n_clusters == 0 && isempty(smld_out.emitters)
    end
    end  # SMLM_TEST_FULL

    @testset "use_3d on 2D data errors" begin
        smld = _make_prec_smld([(0.0, 0.0, 1), (1.0, 1.0, 1)])
        @test_throws ErrorException cluster(smld, PrecisionDBSCANConfig(nsigma = 5.0, use_3d = true))
    end

end
