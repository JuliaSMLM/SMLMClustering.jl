using SMLMClustering
using SMLMData
using Test
using Random

# Helper: build a 2D BasicSMLD of Emitter2DFit{Float64} from (x, y, dataset)
# tuples (microns). frame=1 throughout.
function _make_2d_smld_hdbscan(points::Vector{Tuple{Float64,Float64,Int}};
                              n_datasets::Int = maximum(p[3] for p in points))
    cam = IdealCamera(1:64, 1:64, 0.1)
    emitters = [Emitter2DFit{Float64}(
        x, y, 1000.0, 10.0, 0.01, 0.01, 50.0, 2.0;
        frame = 1, dataset = ds,
    ) for (x, y, ds) in points]
    BasicSMLD(emitters, cam, 1, n_datasets, Dict{String,Any}())
end

function _blob_hdbscan(rng, cx, cy, σ, n; dataset = 1)
    [(cx + σ * randn(rng), cy + σ * randn(rng), dataset) for _ in 1:n]
end

@testset "HDBSCAN backend" begin

    @testset "config construction" begin
        cfg = HDBSCANConfig()
        @test cfg isa AbstractClusterConfig
        @test cfg.min_points == 5
        @test cfg.min_cluster_size === nothing
        @test cfg.knn_graph_k == 30
        @test cfg.cluster_selection_method === :eom
        @test cfg.allow_single_cluster === false
        @test cfg.halo_trim_frac == 0.10
        @test cfg.use_3d === false
        @test cfg.per_dataset === true
        @test cfg.remove_unclustered === false

        cfg2 = HDBSCANConfig(min_points = 10, min_cluster_size = 20,
                             knn_graph_k = 50, cluster_selection_method = :leaf,
                             allow_single_cluster = true,
                             use_3d = true, per_dataset = false,
                             remove_unclustered = true)
        @test cfg2.min_points == 10
        @test cfg2.min_cluster_size == 20
        @test cfg2.knn_graph_k == 50
        @test cfg2.cluster_selection_method === :leaf
        @test cfg2.allow_single_cluster === true
        @test cfg2.use_3d === true
        @test cfg2.per_dataset === false
        @test cfg2.remove_unclustered === true
    end

    @testset "argument validation" begin
        smld = _make_2d_smld_hdbscan([(0.0, 0.0, 1)])
        @test_throws ArgumentError cluster(smld, HDBSCANConfig(min_points = 0))
        @test_throws ArgumentError cluster(smld, HDBSCANConfig(knn_graph_k = 0))
        @test_throws ArgumentError cluster(smld, HDBSCANConfig(min_cluster_size = 1))
        @test_throws ArgumentError cluster(smld, HDBSCANConfig(cluster_selection_method = :bogus))
    end

    @testset "halo trimming keeps a diffuse cluster within its physical extent" begin
        # Regression for the HDBSCAN over-assignment bug: raw HDBSCAN* labels every
        # point under a selected cluster's branch, so a diffuse cluster near the noise
        # floor swept in background out to ~4µm. halo_trim_frac (default 0.10) must keep
        # its members within the cluster's physical extent (~0.7µm here).
        rng = Xoshiro(20260627)
        pts = Tuple{Float64,Float64,Int}[]
        append!(pts, _blob_hdbscan(rng, 1.5, 2.5, 0.03, 400))   # dense
        append!(pts, _blob_hdbscan(rng, 3.5, 2.5, 0.20, 400))   # diffuse, near noise floor
        for _ in 1:150
            push!(pts, (5.0 * rand(rng), 5.0 * rand(rng), 1))   # uniform background
        end
        smld = _make_2d_smld_hdbscan(pts; n_datasets = 1)

        # (label, max radius from the diffuse core) of the cluster covering (3.5, 2.5)
        diffuse_radius = function (o)
            lab = [e.id for e in o.emitters]
            X = [e.x for e in o.emitters]; Y = [e.y for e in o.emitters]
            dl = lab[argmin((X .- 3.5) .^ 2 .+ (Y .- 2.5) .^ 2)]
            mem = findall(==(dl), lab)
            (dl, maximum(sqrt.((X[mem] .- 3.5) .^ 2 .+ (Y[mem] .- 2.5) .^ 2)))
        end

        out, info = cluster(smld, HDBSCANConfig(min_points = 10, min_cluster_size = 25,
                                                per_dataset = false))
        @test info.n_clusters == 2
        dl, rmax = diffuse_radius(out)
        @test dl > 0
        @test rmax < 1.5                     # physical edge ≈ 0.7µm; the bug reached ~4µm

        # halo_trim_frac = 0 disables trimming → the wide halo returns
        out0, info0 = cluster(smld, HDBSCANConfig(min_points = 10, min_cluster_size = 25,
                                                  per_dataset = false, halo_trim_frac = 0.0))
        @test info0.n_clusters == 2
        dl0, rmax0 = diffuse_radius(out0)
        @test dl0 > 0
        @test rmax0 > 2.0

        # out-of-range halo_trim_frac is rejected
        @test_throws ArgumentError cluster(smld,
            HDBSCANConfig(halo_trim_frac = 1.0, per_dataset = false))
    end

end
