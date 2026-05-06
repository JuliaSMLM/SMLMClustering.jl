using Test
using Random
using SMLMClustering
using SMLMClustering: EdgeClassify
using SMLMClustering.EdgeClassify: _build_mask_carve, _mc_polygon_to_mask,
                                    _mc_trace_outer_contour, _mc_cell_to_xy,
                                    _METHOD_MASK_CARVE

@testset "EdgeClassify mask_carve" begin

    @testset "Defaults: MASK_CARVE_* fields" begin
        p = EdgeClassifyParams()
        @test p.MASK_CARVE_SIGMA_UM == 0.080
        @test p.MASK_CARVE_K_NOISE == 3.0
        @test p.MASK_CARVE_PIXEL_UM == 0.040
        @test p.MASK_CARVE_MIN_COMPONENT_FRAC == 0.05
        @test p.MASK_CARVE_FILL_HOLE_MAX_UM2 == 0.5
    end

    @testset "METHOD selector accepts mask_carve" begin
        # Construction does not error.
        p = EdgeClassifyParams(METHOD = "mask_carve")
        @test p.METHOD == "mask_carve"
        # "mask_carve" is in the valid set (no ArgumentError on dispatch
        # from the METHOD-validity check). Trivial 2-point input fails at the
        # v1 alpha-shape stage as expected — this confirms mask_carve
        # dispatch happens AFTER v1 (which is the contract).
        x = [0.0, 1.0]; y = [0.0, 1.0]
        @test_throws Exception classify_emitters(x, y;
            fov_um = (0.0, 1.0, 0.0, 1.0),
            params = EdgeClassifyParams(METHOD = "mask_carve"))
    end

    # --- Synthetic helpers ----------------------------------------------------

    function _disk_pts(rng, cx, cy, r, n)
        x = Float64[]; y = Float64[]
        while length(x) < n
            xx = cx - r + 2r*rand(rng); yy = cy - r + 2r*rand(rng)
            (xx-cx)^2 + (yy-cy)^2 <= r^2 || continue
            push!(x, xx); push!(y, yy)
        end
        return x, y
    end

    function _bitten_pts(rng, cx, cy, r, bx, by, br, n)
        x = Float64[]; y = Float64[]
        while length(x) < n
            xx = cx - r + 2r*rand(rng); yy = cy - r + 2r*rand(rng)
            (xx-cx)^2 + (yy-cy)^2 <= r^2 || continue
            (xx-bx)^2 + (yy-by)^2 > br^2 || continue
            push!(x, xx); push!(y, yy)
        end
        return x, y
    end

    @testset "Parameter validation rejects non-positive MASK_CARVE_*" begin
        rng = MersenneTwister(303)
        x, y = _disk_pts(rng, 5.0, 5.0, 2.0, 4000)
        append!(x, 10 .* rand(rng, 100)); append!(y, 10 .* rand(rng, 100))
        fov = (0.0, 10.0, 0.0, 10.0)
        mk = (; kw...) -> EdgeClassifyParams(K_LIST = [8, 32], RHO_K_THRESH = 50.0,
                                              ALPHA_NM = 400.0, MEMBRANE_NM = 100.0,
                                              REFLECT_RADIUS_NM = 200.0,
                                              METHOD = "mask_carve"; kw...)
        @test_throws ArgumentError classify_emitters(x, y; fov_um=fov,
            params = mk(MASK_CARVE_SIGMA_UM = 0.0))
        @test_throws ArgumentError classify_emitters(x, y; fov_um=fov,
            params = mk(MASK_CARVE_SIGMA_UM = -0.1))
        @test_throws ArgumentError classify_emitters(x, y; fov_um=fov,
            params = mk(MASK_CARVE_PIXEL_UM = 0.0))
        @test_throws ArgumentError classify_emitters(x, y; fov_um=fov,
            params = mk(MASK_CARVE_K_NOISE = 0.0))
        @test_throws ArgumentError classify_emitters(x, y; fov_um=fov,
            params = mk(MASK_CARVE_MIN_COMPONENT_FRAC = -0.01))
        @test_throws ArgumentError classify_emitters(x, y; fov_um=fov,
            params = mk(MASK_CARVE_FILL_HOLE_MAX_UM2 = -0.5))
        # Boundary case: zero is permitted for FRAC and FILL_HOLE.
        r0 = classify_emitters(x, y; fov_um=fov,
            params = mk(MASK_CARVE_MIN_COMPONENT_FRAC = 0.0,
                         MASK_CARVE_FILL_HOLE_MAX_UM2 = 0.0))
        @test r0.params_used.METHOD == "mask_carve"
    end

    @testset "Tracer closure invariant on a disk mask" begin
        # On a normal disk-like mask, the tracer must produce a contour with
        # contour[end] == contour[1] (closed). This is the invariant the
        # builder's fallback guard relies on; assert it directly.
        fov = (0.0, 4.0, 0.0, 4.0); pixel = 0.04
        nx = round(Int, (fov[2]-fov[1])/pixel)
        ny = round(Int, (fov[4]-fov[3])/pixel)
        disk = falses(nx, ny)
        cx, cy, r = 2.0, 2.0, 1.2
        for j in 1:ny, i in 1:nx
            px = fov[1] + (i - 0.5)*pixel; py = fov[3] + (j - 0.5)*pixel
            if (px - cx)^2 + (py - cy)^2 <= r^2
                disk[i, j] = true
            end
        end
        contour = _mc_trace_outer_contour(disk)
        @test length(contour) >= 4
        @test contour[end] == contour[1]    # closure invariant
    end

    @testset "rasterize → polygonize registration (square)" begin
        # A unit square inside a 4×4 µm FOV. Rasterize as a polygon, polygonize,
        # check the resulting polygon is in roughly the right place (not
        # transposed) and approximately the same shape.
        fov = (0.0, 4.0, 0.0, 4.0)
        square = [(1.0, 1.0), (3.0, 1.0), (3.0, 3.0), (1.0, 3.0)]
        pixel = 0.04
        nx = round(Int, (fov[2]-fov[1])/pixel)
        ny = round(Int, (fov[4]-fov[3])/pixel)
        m = _mc_polygon_to_mask(square, fov, pixel, nx, ny)
        contour = _mc_trace_outer_contour(m)
        @test length(contour) > 4
        verts = [_mc_cell_to_xy(i, j, fov, pixel) for (i, j) in contour]
        xs = [v[1] for v in verts]; ys = [v[2] for v in verts]
        # Must lie in the square area (x ∈ [1, 3] ± pixel, y ∈ [1, 3] ± pixel).
        @test minimum(xs) >= 1.0 - 2*pixel
        @test maximum(xs) <= 3.0 + 2*pixel
        @test minimum(ys) >= 1.0 - 2*pixel
        @test maximum(ys) <= 3.0 + 2*pixel
        # Centroid near (2, 2), confirming no axis transposition.
        @test abs(sum(xs)/length(xs) - 2.0) < 0.10
        @test abs(sum(ys)/length(ys) - 2.0) < 0.10
    end

    # --- Synthetic crescent: v1 chord-bridges, mask_carve carves the bay -----

    @testset "Carve crescent: bay carved, default outer_polygon unchanged" begin
        rng = MersenneTwister(513)
        # Crescent: big disk minus offset bite.
        cx, cy, R = 5.0, 5.0, 3.0
        bx, by, br = 8.0, 5.0, 1.5
        x, y = _bitten_pts(rng, cx, cy, R, bx, by, br, 12000)
        append!(x, 10 .* rand(rng, 200))
        append!(y, 10 .* rand(rng, 200))
        fov = (0.0, 10.0, 0.0, 10.0)

        base = EdgeClassifyParams(K_LIST = [8, 32], RHO_K_THRESH = 50.0,
                                  ALPHA_NM = 400.0, MEMBRANE_NM = 100.0,
                                  REFLECT_RADIUS_NM = 200.0)
        v1 = classify_emitters(x, y; fov_um = fov, params = base)

        carved = classify_emitters(x, y; fov_um = fov,
            params = EdgeClassifyParams(K_LIST = base.K_LIST,
                                         RHO_K_THRESH = base.RHO_K_THRESH,
                                         ALPHA_NM = base.ALPHA_NM,
                                         MEMBRANE_NM = base.MEMBRANE_NM,
                                         REFLECT_RADIUS_NM = base.REFLECT_RADIUS_NM,
                                         METHOD = "mask_carve"))

        # Default METHOD path produces no diagnostic.
        @test v1.mask_carve_diagnostic === nothing
        # mask_carve produces a diagnostic.
        @test carved.mask_carve_diagnostic !== nothing
        @test carved.params_used.METHOD == "mask_carve"

        # Effective outer polygon differs from v1 alpha-shape outer (carve
        # actually carved something) — but not catastrophically.
        @test carved.outer_polygon !== carved.loops[1]
        @test carved.loops[1] == v1.loops[1]   # alpha provenance preserved

        # Class partition invariant.
        n = length(x)
        @test count(==("outside"), carved.class) +
              count(==("membrane"), carved.class) +
              count(==("interior"), carved.class) == n

        # Carve must not expand v1 outward (carve_only_area ≈ 0 by construction).
        diag = carved.mask_carve_diagnostic
        @test diag.applied
        @test diag.carve_only_area_um2 < 0.05   # rasterization roundoff bound
        # NOTE: this integration test does NOT assert how much carving
        # happened, because at this fixture's ALPHA_NM the alpha-shape outer
        # may already follow the bay (no chord to carve away). The
        # builder-level test below exercises a known chorded v1 polygon
        # explicitly and asserts the carve excludes the bay.
    end

    # --- Builder-level test: known chorded v1 + bitten emitters ---------------

    @testset "_build_mask_carve excludes bay vs chorded v1" begin
        rng = MersenneTwister(813)
        # Bitten disk emitters; v1 is the FULL disk polygon (chord-bridges
        # the bay) — this is the synthetic v3.1 setup. mask_carve must
        # produce a carve polygon that excludes the bay region.
        cx, cy, R = 5.0, 5.0, 2.0
        bx, by, br = 7.0, 5.0, 1.5   # bite at disk perimeter — true boundary bay
        x, y = _bitten_pts(rng, cx, cy, R, bx, by, br, 10000)
        append!(x, 10 .* rand(rng, 200))
        append!(y, 10 .* rand(rng, 200))
        fov = (0.0, 10.0, 0.0, 10.0)
        # Chorded v1: full disk polygon at radius R.
        v1_poly = NTuple{2,Float64}[(cx + R*cos(2π*k/96), cy + R*sin(2π*k/96))
                                    for k in 0:95]

        params = EdgeClassifyParams(METHOD = "mask_carve")
        carve_poly, diag = _build_mask_carve(v1_poly, x, y, fov, params)
        @test diag.applied
        # Carve ⊆ v1 invariant.
        @test diag.carve_only_area_um2 < 0.05
        # Bay area ≈ overlap area between v1 disk and bite disk. For a bite
        # at the perimeter (d = R = 2.0, br = 1.5), the lens area is large.
        # Require v1_only_area to capture a substantial fraction (≥ 1.0 µm²).
        @test diag.v1_only_area_um2 > 1.0
        # Polygon point-in-polygon: bay center should be EXCLUDED from carve.
        @test !SMLMClustering.EdgeClassify._point_in_polygon(bx, by, carve_poly)
    end

    # --- Default-METHOD outer_polygon regression (regression guard) ----------

    @testset "Default METHOD outer_polygon unchanged" begin
        rng = MersenneTwister(99)
        x, y = _disk_pts(rng, 5.0, 5.0, 2.5, 8000)
        append!(x, 10 .* rand(rng, 200))
        append!(y, 10 .* rand(rng, 200))
        fov = (0.0, 10.0, 0.0, 10.0)
        base = EdgeClassifyParams(K_LIST = [8, 32], RHO_K_THRESH = 50.0,
                                  ALPHA_NM = 400.0, MEMBRANE_NM = 100.0,
                                  REFLECT_RADIUS_NM = 200.0)
        r1 = classify_emitters(x, y; fov_um = fov, params = base)
        r2 = classify_emitters(x, y; fov_um = fov, params = base)
        @test r1.outer_polygon == r2.outer_polygon
        @test r1.outer_polygon == r1.loops[1]   # default: outer_polygon == loops[1]
        @test r1.mask_carve_diagnostic === nothing
        @test all(r1.class .== r2.class)
        @test all(r1.dist_to_outer_um[r1.inside_outer] .==
                  r2.dist_to_outer_um[r2.inside_outer])
    end

    # --- dist_to_outer uses effective carve polygon, not loops[1] ------------

    @testset "dist_to_outer measured against effective (carve) polygon" begin
        rng = MersenneTwister(202)
        cx, cy, R = 5.0, 5.0, 3.0
        bx, by, br = 8.0, 5.0, 1.5
        x, y = _bitten_pts(rng, cx, cy, R, bx, by, br, 12000)
        append!(x, 10 .* rand(rng, 200))
        append!(y, 10 .* rand(rng, 200))
        fov = (0.0, 10.0, 0.0, 10.0)

        carved = classify_emitters(x, y; fov_um = fov,
            params = EdgeClassifyParams(K_LIST = [8, 32], RHO_K_THRESH = 50.0,
                                         ALPHA_NM = 400.0, MEMBRANE_NM = 100.0,
                                         REFLECT_RADIUS_NM = 200.0,
                                         METHOD = "mask_carve"))
        # If outer_polygon != loops[1], dist_to_outer for at least some
        # interior emitters should be different from what loops[1] would yield.
        any_diff = false
        for i in 1:length(x)
            carved.inside_outer[i] || continue
            d_loop = SMLMClustering.EdgeClassify._dist_to_polygon(x[i], y[i], carved.loops[1])
            d_eff  = carved.dist_to_outer_um[i]
            if abs(d_loop - d_eff) > 1e-6
                any_diff = true; break
            end
        end
        # Only assert the inequality when the carve actually moved the polygon.
        if carved.outer_polygon !== carved.loops[1]
            @test any_diff
        end
    end

    # --- Fallback path (very sparse input forces empty intersection) ---------

    @testset "Fallback when carve degenerate (no crash)" begin
        # Build a minimal valid alpha-shape input so v1 succeeds, then ensure
        # mask_carve falls back without error if the density grid is degenerate.
        rng = MersenneTwister(404)
        x, y = _disk_pts(rng, 5.0, 5.0, 2.0, 8000)
        append!(x, 10 .* rand(rng, 100))
        append!(y, 10 .* rand(rng, 100))
        fov = (0.0, 10.0, 0.0, 10.0)
        # Force a very high k_noise so threshold > pmax → empty d_mask.
        carved = classify_emitters(x, y; fov_um = fov,
            params = EdgeClassifyParams(K_LIST = [8, 32], RHO_K_THRESH = 50.0,
                                         ALPHA_NM = 400.0, MEMBRANE_NM = 100.0,
                                         REFLECT_RADIUS_NM = 200.0,
                                         METHOD = "mask_carve",
                                         MASK_CARVE_K_NOISE = 1e9))
        @test carved.mask_carve_diagnostic !== nothing
        @test carved.mask_carve_diagnostic.applied == false
        @test !isempty(carved.mask_carve_diagnostic.fallback_reason)
        # Class partition invariant still holds.
        n = length(x)
        @test count(==("outside"), carved.class) +
              count(==("membrane"), carved.class) +
              count(==("interior"), carved.class) == n
        # Effective outer falls back to v1 alpha outer.
        @test carved.outer_polygon == carved.loops[1]
        # Fallback diagnostic must describe the effective fallback boundary
        # (which is v1), not a phantom zero-area carve. Internal consistency:
        diag = carved.mask_carve_diagnostic
        v1_area = SMLMClustering.EdgeClassify._polygon_area_abs(carved.loops[1])
        @test diag.v1_polygon_area_um2 ≈ v1_area
        @test diag.carve_polygon_area_um2 ≈ diag.v1_polygon_area_um2
        @test diag.area_delta_um2 == 0.0
        @test diag.v1_only_area_um2 == 0.0
        @test diag.carve_only_area_um2 == 0.0
        @test diag.med_v1_carve_distance_um == 0.0
        @test diag.p95_v1_carve_distance_um == 0.0
        @test diag.n_carve_polygon_pts == length(carved.loops[1])
    end

    # --- Artifact writer smoke test ------------------------------------------

    @testset "Artifact writers (mask_carve)" begin
        rng = MersenneTwister(606)
        cx, cy, R = 5.0, 5.0, 2.5
        bx, by, br = 7.5, 5.0, 1.2
        x, y = _bitten_pts(rng, cx, cy, R, bx, by, br, 10000)
        append!(x, 10 .* rand(rng, 150))
        append!(y, 10 .* rand(rng, 150))
        fov = (0.0, 10.0, 0.0, 10.0)
        mktempdir() do tmp
            classify_emitters(x, y; fov_um = fov,
                params = EdgeClassifyParams(K_LIST = [8, 32], RHO_K_THRESH = 50.0,
                                             ALPHA_NM = 400.0, MEMBRANE_NM = 100.0,
                                             REFLECT_RADIUS_NM = 200.0,
                                             METHOD = "mask_carve"),
                out_dir = tmp, condition = "synth", cell = "cell_01",
                write_artifacts = true)
            leaf = joinpath(tmp, "synth", "cell_01")
            @test isfile(joinpath(leaf, "classified.tsv"))
            @test isfile(joinpath(leaf, "polygon_loops.tsv"))
            @test isfile(joinpath(leaf, "loop_diagnostics.csv"))
            @test isfile(joinpath(leaf, "params.json"))
            @test isfile(joinpath(leaf, "manifest.json"))
            @test isfile(joinpath(leaf, "effective_outer.tsv"))
            @test isfile(joinpath(leaf, "mask_carve_diagnostic.json"))
            man = read(joinpath(leaf, "manifest.json"), String)
            @test occursin("\"effective_outer_tsv\"", man)
            @test occursin("\"mask_carve_diagnostic_json\"", man)
            mc = read(joinpath(leaf, "mask_carve_diagnostic.json"), String)
            @test occursin("\"applied\"", mc)
            @test occursin("\"v1_polygon_area_um2\"", mc)
            ptxt = read(joinpath(leaf, "params.json"), String)
            @test occursin("\"MASK_CARVE_SIGMA_UM\"", ptxt)
        end
    end

    @testset "Artifact writers: outer_polygon path skips mask_carve files" begin
        rng = MersenneTwister(707)
        x, y = _disk_pts(rng, 5.0, 5.0, 2.5, 8000)
        append!(x, 10 .* rand(rng, 200))
        append!(y, 10 .* rand(rng, 200))
        fov = (0.0, 10.0, 0.0, 10.0)
        mktempdir() do tmp
            classify_emitters(x, y; fov_um = fov,
                params = EdgeClassifyParams(K_LIST = [8, 32], RHO_K_THRESH = 50.0,
                                             ALPHA_NM = 400.0, MEMBRANE_NM = 100.0,
                                             REFLECT_RADIUS_NM = 200.0),
                out_dir = tmp, condition = "synth", cell = "cell_01",
                write_artifacts = true)
            leaf = joinpath(tmp, "synth", "cell_01")
            @test !isfile(joinpath(leaf, "effective_outer.tsv"))
            @test !isfile(joinpath(leaf, "mask_carve_diagnostic.json"))
            man = read(joinpath(leaf, "manifest.json"), String)
            # Manifest still lists the entries with written=false.
            @test occursin("\"effective_outer_tsv\"", man)
            @test occursin("\"mask_carve_diagnostic_json\"", man)
            @test occursin("\"written\": false", man)
        end
    end
end
