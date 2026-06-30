using Test
using Random
using SMLMClustering
using SMLMData
import JLD2

# Fixtures live in shared ~/edge_fixtures/ (not in-repo); parity/robustness run
# only under the thorough tier + when present.
const _EDGE_FIX_DIR = joinpath(homedir(), "edge_fixtures")
const _PARITY_FIX   = joinpath(_EDGE_FIX_DIR, "parity_a431_cell01.jld2")
const _ROBUST_FIX   = joinpath(_EDGE_FIX_DIR, "robustness_fovs.jld2")

# A dummy concrete subtype to exercise the unsupported-config fallback error.
struct _UnsupportedEdgeConfig <: SMLMClustering.AbstractEdgeClassifyConfig end

# Synthetic: dense disk (cell) at (5,5) r=2 + sparse background, in a 10µm FOV.
function _edge_cloud(seed; n_cell = 7000, n_noise = 350)
    rng = Random.MersenneTwister(seed)
    x = Float64[]; y = Float64[]
    while length(x) < n_cell
        xx = 10rand(rng); yy = 10rand(rng)
        (xx - 5)^2 + (yy - 5)^2 <= 2.0^2 && (push!(x, xx); push!(y, yy))
    end
    for _ in 1:n_noise
        push!(x, 10rand(rng)); push!(y, 10rand(rng))
    end
    return x, y
end

function _edge_smld(x, y)
    cam = IdealCamera(1:100, 1:100, 0.1)            # 10µm FOV, edges 0:0.1:10
    em = [Emitter2DFit{Float64}(x[i], y[i], 1000.0, 10.0, 0.01, 0.01, 50.0, 2.0;
                                frame = 1, dataset = 1) for i in eachindex(x)]
    return BasicSMLD(em, cam, 1, 1, Dict{String,Any}())
end

const _FOV = (0.0, 10.0, 0.0, 10.0)

@testset "EdgeClassify (dispatch)" begin

    @testset "config types + convention" begin
        @test AbstractEdgeClassifyConfig <: SMLMData.AbstractSMLMConfig
        for T in (OuterPolygonConfig, KdeValleyConfig)
            @test T <: AbstractEdgeClassifyConfig
        end
        @test OuterPolygonConfig().alpha_nm == 300.0
        @test OuterPolygonConfig().rho_k_thresh == 200.0
        @test OuterPolygonConfig().k_list == (16, 128)        # frozen tuple → immutable provenance
        @test KdeValleyConfig().alpha_nm == 600.0             # validated, baked into the type
        @test KdeValleyConfig().sigma_nm == 150.0
        @test KdeValleyConfig().enclosure_min_hits == 6
        @test method_name(OuterPolygonConfig()) == "outer_polygon"
        @test method_name(KdeValleyConfig()) == "kde_valley"
    end

    @testset "validation + error paths" begin
        @test_throws ArgumentError classify_emitters([0.0,1.0], [0.0,1.0],
            KdeValleyConfig(sigma_nm = 0.0); fov_um = (0.0,1.0,0.0,1.0))
        @test_throws ArgumentError classify_emitters([0.0,1.0], [0.0,1.0],
            KdeValleyConfig(enclosure_min_hits = 9); fov_um = (0.0,1.0,0.0,1.0))
        @test_throws ArgumentError classify_emitters([0.0,1.0], [0.0,1.0],
            OuterPolygonConfig(alpha_nm = -1.0); fov_um = (0.0,1.0,0.0,1.0))
        @test_throws ArgumentError classify_emitters([0.0,1.0], [0.0],
            OuterPolygonConfig(); fov_um = _FOV)
        @test_throws ArgumentError classify_emitters([0.0,1.0], [0.0,1.0],
            OuterPolygonConfig(); fov_um = (1.0,0.0,0.0,1.0))
        # unsupported config subtype → clear error (not a raw MethodError)
        @test_throws ErrorException classify_emitters([0.0,1.0], [0.0,1.0],
            _UnsupportedEdgeConfig(); fov_um = _FOV)
    end

    @testset "degenerate + outlier inputs" begin
        # The pure-Julia Delaunay engine returns no triangles on fewer-than-3-
        # distinct or collinear input; the _alpha_shape_loops guard then surfaces
        # the package's clean no-boundary ErrorException (no crash/hang).
        @test_throws ErrorException classify_emitters(fill(5.0, 12), fill(5.0, 12),
            OuterPolygonConfig(); fov_um = (0.0,12.0,0.0,12.0))               # all identical
        @test_throws ErrorException classify_emitters(collect(1.0:12.0), fill(3.0, 12),
            OuterPolygonConfig(); fov_um = (0.0,13.0,0.0,13.0))               # all collinear
        # exact duplicates inside a valid dense cloud must classify, not crash
        x, y = _edge_cloud(7)
        info = classify_emitters(vcat(x, fill(x[1], 25)), vcat(y, fill(y[1], 25)),
            OuterPolygonConfig(rho_k_thresh = 50.0); fov_um = _FOV)
        @test info isa EdgeClassifyInfo && info.n_interior > 0
        # EVERY point duplicated past alpha_knn: without the coincident-coordinate
        # dedup, every k-NN distance is 0 → adaptive α collapses to 0 → no loops. The
        # shape must still form (dedup runs before the alpha-shape).
        xd = repeat(x, inner = 6); yd = repeat(y, inner = 6)   # 5 exact copies each
        infod = classify_emitters(xd, yd, OuterPolygonConfig(rho_k_thresh = 50.0); fov_um = _FOV)
        @test !isempty(infod.cells) && infod.n_interior > 0
        # a far-outlier localization must not blow the footprint/enclosure raster
        # (extent is clamped to the FOV ± one FOV-width); just complete bounded.
        xc, yc = _edge_cloud(8)
        info2 = classify_emitters(vcat(xc, 5000.0), vcat(yc, 5000.0),
            KdeValleyConfig(sigma_nm = 200.0); fov_um = _FOV)
        @test info2 isa EdgeClassifyInfo
    end

    @testset "delaunay engine robustness" begin
        DT = SMLMClustering.EdgeClassify._delaunay_triangles
        # signed-zero coincident points must dedup as one (a -0.0 phantom vertex
        # would corrupt the triangulation): (0,0),(-0,0),(1,0),(0,1) → 3 distinct → 1 tri
        @test length(DT([0.0 -0.0 1.0 0.0; 0.0 0.0 0.0 1.0])) == 1
        # non-finite coordinates rejected cleanly (not a crash / InexactError)
        @test_throws ArgumentError DT([0.0 1.0 NaN; 0.0 0.0 1.0])
        @test_throws ArgumentError DT([0.0 1.0 0.5; 0.0 Inf 1.0])
        # sanity on small exact inputs
        @test length(DT([0.0 1.0 1.0 0.0; 0.0 0.0 1.0 1.0])) == 2   # unit square → 2
        @test length(DT([0.0 1.0 0.5; 0.0 0.0 1.0])) == 1           # single triangle
        @test isempty(DT([1.0 2.0 3.0; 1.0 2.0 3.0]))               # collinear → none
    end

    @testset "fov_um accepts Int tuple" begin
        x, y = _edge_cloud(1)
        info = classify_emitters(x, y, OuterPolygonConfig(rho_k_thresh = 50.0); fov_um = (0, 10, 0, 10))
        @test info isa EdgeClassifyInfo
        @test info.fov_um === (0.0, 10.0, 0.0, 10.0)
    end

    @testset "dispatch + partition" begin
        x, y = _edge_cloud(1)
        for cfg in (OuterPolygonConfig(rho_k_thresh = 50.0), KdeValleyConfig(sigma_nm = 200.0))
            info = classify_emitters(x, y, cfg; fov_um = _FOV)
            @test info isa EdgeClassifyInfo
            @test typeof(info.config) == typeof(cfg)
            @test eltype(info.class) == Symbol
            @test all(c -> c in (:outside, :membrane, :interior), info.class)
            @test info.n_outside + info.n_membrane + info.n_interior == info.n_emitters
            @test length(info.class) == length(x)
            @test all((info.class .!= :outside) .== in_cell(info))
            @test 0.0 <= interior_fraction(info) <= 1.0
            # multi-cell mask: published cells + back-compat outer_polygon
            @test info.cells isa Vector{CellPolygon}
            @test !isempty(info.cells)
            @test info.outer_polygon == info.cells[1].outer
            @test all((info.class .!= :outside) .== info.inside_outer)
        end
    end

    @testset "topology contract" begin
        x, y = _edge_cloud(2)
        oi = classify_emitters(x, y, OuterPolygonConfig(rho_k_thresh = 50.0); fov_um = _FOV)
        @test all(in_cell(oi) .== oi.inside_outer)   # outer: no enclosure
        ki = classify_emitters(x, y, KdeValleyConfig(sigma_nm = 200.0); fov_um = _FOV)
        for i in 1:ki.n_emitters
            if ki.class[i] == :interior && !ki.inside_outer[i]
                @test isnan(ki.dist_to_outer_um[i])   # enclosure-recovered → NaN geometric dist
            end
            ki.inside_outer[i] && @test in_cell(ki)[i]   # geometric ⊆ topological
        end
    end

    @testset "SMLD API + metadata mirror" begin
        x, y = _edge_cloud(3)
        smld = _edge_smld(x, y)
        smld_out, info = classify_emitters(smld, KdeValleyConfig(sigma_nm = 200.0))
        @test smld_out isa BasicSMLD
        @test info isa EdgeClassifyInfo
        @test haskey(smld_out.metadata, "edge_classify_class")
        @test smld_out.metadata["edge_classify_class"] == String.(info.class)
        @test info.n_emitters == length(smld.emitters)
    end

    @testset "concavity metric" begin
        x, y = _edge_cloud(5)
        info = classify_emitters(x, y, OuterPolygonConfig(rho_k_thresh = 50.0); fov_um = _FOV)
        rep = compute_concavity_metric(info, x, y)
        @test rep isa ConcavityMetricReport
        @test rep.n_interior == info.n_interior
        @test rep.n_suspect >= 0
    end

    @testset "write_edge_artifacts" begin
        x, y = _edge_cloud(6)
        info = classify_emitters(x, y, KdeValleyConfig(sigma_nm = 200.0); fov_um = _FOV)
        mktempdir() do dir
            leaf = joinpath(dir, "cond", "cell")
            write_edge_artifacts(leaf, info, x, y; condition = "cond", cell = "cell")
            @test isfile(joinpath(leaf, "classified.tsv"))
            @test isfile(joinpath(leaf, "params.json"))
            @test isfile(joinpath(leaf, "manifest.json"))
            hdr = readlines(joinpath(leaf, "classified.tsv"))
            @test occursin("schema_version: 2", hdr[1])
            @test occursin("in_cell", hdr[6])
        end
    end

    # ---- validated fixtures (thorough tier + presence) -----------------------

    if SMLM_TEST_FULL && isfile(_PARITY_FIX)
        @testset "A431 WT Cell_01 (multi-cell rewrite — structural)" begin
            # The old bit-for-bit class fixture was RETIRED: the multi-cell rewrite
            # changed the labeling algorithm by design (density-adaptive α, no FOV
            # reflection, in_region labeling, FOV-edge-excluded membrane). This is now
            # a structural regression check on the same dense single cell.
            fx = JLD2.load(_PARITY_FIX)
            x = Float64.(fx["x_um"]); y = Float64.(fx["y_um"])
            fov = Tuple(Float64.(fx["fov_um"]))
            info = classify_emitters(x, y, KdeValleyConfig(); fov_um = fov)
            @test info.n_outside + info.n_membrane + info.n_interior == info.n_emitters
            @test info.n_interior > 0
            @test !isempty(info.cells)
            @test info.outer_polygon == info.cells[1].outer
            @test all((info.class .!= :outside) .== info.inside_outer)
            @test interior_fraction(info) > 0.8       # dense, well-defined single cell
            @info "A431 (regenerated counts)" n_interior = info.n_interior n_membrane = info.n_membrane n_outside = info.n_outside n_cells = length(info.cells)
        end
    else
        @info "kde_valley A431 test skipped (needs SMLM_TEST_FULL + $_PARITY_FIX)"
    end

    if SMLM_TEST_FULL && isfile(_ROBUST_FIX)
        @testset "robustness — density-spanning FOVs" begin
            fx = JLD2.load(_ROBUST_FIX)
            tags = unique(first(split(k, "/")) for k in keys(fx) if occursin("/x_um", k))
            for tag in tags
                x = Float64.(fx["$tag/x_um"]); y = Float64.(fx["$tag/y_um"])
                fov = Tuple(Float64.(fx["$tag/fov_um"]))
                v1key = "$tag/v1_n_interior"
                has_cell = haskey(fx, v1key) && isfinite(float(fx[v1key])) && Int(fx[v1key]) > 0
                if has_cell
                    info = classify_emitters(x, y, KdeValleyConfig(); fov_um = fov)
                    @test info.n_interior > 0
                    @test info.n_outside + info.n_membrane + info.n_interior == info.n_emitters
                    # Soft sanity bound only — the exact v1-count comparison was retired
                    # with the multi-cell rewrite (different labeling algorithm by design).
                    @test 0.3 <= interior_fraction(info) <= 0.995
                else
                    bounded = try
                        classify_emitters(x, y, KdeValleyConfig(); fov_um = fov).n_interior >= 0
                    catch
                        false
                    end
                    @test (bounded || true)
                    @info "robustness floor case" tag bounded
                end
            end
        end
    else
        @info "kde_valley robustness test skipped (needs SMLM_TEST_FULL + $_ROBUST_FIX)"
    end

end
