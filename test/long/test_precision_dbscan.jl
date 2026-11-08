using SMLMClustering
# The reuse primitive is public but not exported — bring the names into scope explicitly.
using SMLMClustering: build_precision_neighbor_graph, precision_dbscan_labels,
                      precision_dbscan_labels!, PrecisionNeighborGraph
using SMLMData
using Clustering        # to cross-check min_points against classical DBSCAN minPts
using Test
using Random

# Same 2D SMLD builder shape as `_make_2d_smld` in long/utils/helpers.jl (precision σ_x = σ_y = 0.01 μm).
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

@testset "Precision-weighted DBSCAN (thorough)" begin

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

end
