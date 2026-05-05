using Test
using Random
using SMLMClustering
using SMLMClustering: EdgeClassify

@testset "EdgeClassify v1" begin

    @testset "EdgeClassifyParams defaults" begin
        p = EdgeClassifyParams()
        @test p.K_LIST == [16, 128]
        @test p.RHO_K_THRESH == 200.0
        @test p.ALPHA_NM == 300.0
        @test p.REFLECT_RADIUS_NM == 1500.0
        @test p.MEMBRANE_NM == 100.0
        @test p.FOV_TRUNC_TOL_NM == 150.0
        @test p.METHOD == "outer_polygon"
        @test p.GRID_PX_NM == 50.0
        @test p.GRID_SMOOTH_NM == 80.0
        @test p.GRID_MASK_Q == 0.03
        @test p.GRID_MASK_PEAK_FRAC == 0.26
        @test p.GRID_OUTER_BUFFER_NM == 800.0
        @test p.CONCAVITY_METRIC_BUFFER_NM == 2000.0
    end

    @testset "METHOD selector" begin
        x = [0.0, 1.0]; y = [0.0, 1.0]
        # Unknown method errors
        @test_throws ArgumentError classify_emitters(x, y;
            fov_um = (0.0, 1.0, 0.0, 1.0),
            params = EdgeClassifyParams(METHOD = "bogus"))
        # concave_refined reserved → not implemented yet on baseline path
        @test_throws ArgumentError classify_emitters(x, y;
            fov_um = (0.0, 1.0, 0.0, 1.0),
            params = EdgeClassifyParams(METHOD = "concave_refined"))
    end

    @testset "grid_hybrid synthetic GT" begin
        function membrane_prf(pred, truth)
            tp = count(i -> pred[i] == "membrane" && truth[i] == "membrane", eachindex(pred))
            fp = count(i -> pred[i] == "membrane" && truth[i] != "membrane", eachindex(pred))
            fn = count(i -> pred[i] != "membrane" && truth[i] == "membrane", eachindex(pred))
            precision = tp / max(tp + fp, 1)
            recall = tp / max(tp + fn, 1)
            f1 = (precision + recall) == 0 ? 0.0 : 2 * precision * recall / (precision + recall)
            return precision, recall, f1
        end

        function sample_crescent(seed; n_cell = 16000, n_noise = 800)
            rng = Random.MersenneTwister(seed)
            cx_big, cy_big, R_big = 5.0, 5.0, 3.0
            cx_sm, cy_sm, R_sm = 7.5, 5.0, 2.2
            x = Float64[]; y = Float64[]
            while length(x) < n_cell
                xx = cx_big - R_big + 2 * R_big * rand(rng)
                yy = cy_big - R_big + 2 * R_big * rand(rng)
                in_big = (xx - cx_big)^2 + (yy - cy_big)^2 <= R_big^2
                in_small = (xx - cx_sm)^2 + (yy - cy_sm)^2 <= R_sm^2
                if in_big && !in_small
                    push!(x, xx); push!(y, yy)
                end
            end
            append!(x, 10 .* rand(rng, n_noise))
            append!(y, 10 .* rand(rng, n_noise))
            truth = Vector{String}(undef, length(x))
            for i in eachindex(x)
                in_big = (x[i] - cx_big)^2 + (y[i] - cy_big)^2 <= R_big^2
                in_small = (x[i] - cx_sm)^2 + (y[i] - cy_sm)^2 <= R_sm^2
                if !(in_big && !in_small)
                    truth[i] = "outside"
                else
                    d_big = abs(hypot(x[i] - cx_big, y[i] - cy_big) - R_big)
                    d_bay = abs(hypot(x[i] - cx_sm, y[i] - cy_sm) - R_sm)
                    truth[i] = min(d_big, d_bay) <= 0.10 ? "membrane" : "interior"
                end
            end
            return x, y, truth
        end

        x, y, truth = sample_crescent(72)
        base_params = EdgeClassifyParams(K_LIST = [8, 32], RHO_K_THRESH = 50.0,
                                         ALPHA_NM = 400.0, MEMBRANE_NM = 100.0,
                                         REFLECT_RADIUS_NM = 200.0)
        v1 = classify_emitters(x, y; fov_um = (0.0, 10.0, 0.0, 10.0),
                               params = base_params)
        hybrid = classify_emitters(x, y; fov_um = (0.0, 10.0, 0.0, 10.0),
                                   params = EdgeClassifyParams(
                                       K_LIST = base_params.K_LIST,
                                       RHO_K_THRESH = base_params.RHO_K_THRESH,
                                       ALPHA_NM = base_params.ALPHA_NM,
                                       MEMBRANE_NM = base_params.MEMBRANE_NM,
                                       REFLECT_RADIUS_NM = base_params.REFLECT_RADIUS_NM,
                                       METHOD = "grid_hybrid"))
        vp, vr, vf = membrane_prf(v1.class, truth)
        hp, hr, hf = membrane_prf(hybrid.class, truth)
        @test count(==("outside"), hybrid.class) == count(==("outside"), v1.class)
        @test hr > vr + 0.10
        @test hf > vf + 0.05
        @test hp > 0.75
    end

    @testset "fov_um validation" begin
        x = [0.0]; y = [0.0]
        @test_throws ArgumentError classify_emitters(x, y; fov_um = (1.0, 0.0, 0.0, 1.0))
        @test_throws ArgumentError classify_emitters(x, y; fov_um = (0.0, 1.0, 1.0, 0.0))
        # mismatched lengths
        @test_throws ArgumentError classify_emitters([0.0, 1.0], [0.0]; fov_um = (0.0, 1.0, 0.0, 1.0))
    end

    @testset "synthetic disk: outer-only classifier" begin
        # Synthetic dense disk centered in a 10x10 µm FOV, no FOV truncation.
        rng = Random.MersenneTwister(42)
        cx = cy = 5.0; R = 2.5
        npts = 12000
        # Sample uniformly inside disk: r = R * sqrt(u), θ = 2π v.
        u = rand(rng, npts); v = rand(rng, npts)
        r = R .* sqrt.(u); θ = 2π .* v
        x = cx .+ r .* cos.(θ); y = cy .+ r .* sin.(θ)
        # Add a few sparse "outside" emitters far from disk.
        push!(x, 0.5, 0.5, 9.5); push!(y, 0.5, 9.5, 9.5)

        result = classify_emitters(x, y; fov_um = (0.0, 10.0, 0.0, 10.0),
                                   params = EdgeClassifyParams(
                                       K_LIST = [8, 32], RHO_K_THRESH = 50.0,
                                       ALPHA_NM = 600.0, MEMBRANE_NM = 150.0,
                                       REFLECT_RADIUS_NM = 200.0))

        # Class partition invariant.
        @test length(result.class) == length(x)
        cset = Set(result.class)
        @test cset ⊆ Set(["outside","membrane","interior"])
        @test sum(==("outside"), result.class) +
              sum(==("membrane"), result.class) +
              sum(==("interior"), result.class) == length(x)

        # The 3 hand-placed corner emitters should classify as outside.
        @test all(==("outside"), result.class[end-2:end])

        # Most disk emitters should be interior (membrane is a thin annulus).
        n_disk = npts
        n_interior_disk = sum(==("interior"), result.class[1:n_disk])
        @test n_interior_disk / n_disk > 0.7

        # Outer polygon must be non-trivial.
        @test length(result.outer_polygon) >= 4

        # dist_to_outer_um is NaN exactly where inside_outer is false.
        for i in 1:length(x)
            if result.inside_outer[i]
                @test !isnan(result.dist_to_outer_um[i])
            else
                @test isnan(result.dist_to_outer_um[i])
            end
        end

        # loop_diagnostics: includes loop_id == 1 (outer)
        ldiag = result.loop_diagnostics
        @test !isempty(ldiag)
        @test ldiag[1].loop_id == 1
        @test ldiag[1].heuristic_type == "outer"
        @test 0.0 <= ldiag[1].frac_in_fov <= 1.0
        @test 0.0 <= ldiag[1].frac_dense <= 1.0
    end

    @testset "artifacts + manifest schema" begin
        rng = Random.MersenneTwister(1)
        cx = cy = 5.0; R = 2.5
        npts = 8000
        u = rand(rng, npts); v = rand(rng, npts)
        r = R .* sqrt.(u); θ = 2π .* v
        x = cx .+ r .* cos.(θ); y = cy .+ r .* sin.(θ)

        mktempdir() do tdir
            result = classify_emitters(x, y;
                fov_um = (0.0, 10.0, 0.0, 10.0),
                params = EdgeClassifyParams(K_LIST = [8, 32], RHO_K_THRESH = 50.0,
                                           ALPHA_NM = 600.0, REFLECT_RADIUS_NM = 200.0),
                out_dir = tdir, condition = "TEST", cell = "cell_X",
                write_artifacts = true)

            leaf = joinpath(tdir, "TEST", "cell_X")
            for f in ("classified.tsv","polygon_loops.tsv","loop_diagnostics.csv",
                      "params.json","manifest.json")
                @test isfile(joinpath(leaf, f))
            end

            # classified.tsv: header + one row per emitter
            lines = readlines(joinpath(leaf, "classified.tsv"))
            header_lines = count(l -> startswith(l, "#"), lines)
            colhdr = findfirst(l -> startswith(l, "emitter_id"), lines)
            @test colhdr !== nothing
            data_rows = length(lines) - colhdr
            @test data_rows == result.n_emitters

            # polygon_loops.tsv columns
            ploops = readlines(joinpath(leaf, "polygon_loops.tsv"))
            @test any(l -> startswith(l, "loop_id\tvertex_id\tx_um\ty_um"), ploops)

            # loop_diagnostics.csv schema 2 + used_in_outer before heuristic_type LAST
            lcsv = readlines(joinpath(leaf, "loop_diagnostics.csv"))
            schema_line = first(filter(l -> startswith(l, "# schema_version"), lcsv))
            @test occursin("schema_version: 2", schema_line)
            colhdr = first(filter(l -> startswith(l, "loop_id,"), lcsv))
            cols = split(colhdr, ",")
            @test cols == ["loop_id","vertex_count","area_um2","n_emitters_inside",
                           "frac_in_fov","frac_dense","median_rhoK",
                           "used_in_outer","heuristic_type"]
            @test cols[end] == "heuristic_type"
            @test cols[end-1] == "used_in_outer"

            # manifest.json basic well-formedness
            mtxt = read(joinpath(leaf, "manifest.json"), String)
            @test occursin("\"schema_version\"", mtxt)
            @test occursin("\"classified_tsv\"", mtxt)
            @test occursin("\"manifest.json\"", mtxt) == false   # not self-listed
            @test occursin("\"out_dir\"", mtxt)
            @test occursin("\"leaf_dir\"", mtxt)

            # params.json carries the actual params used (incl. METHOD,
            # CONCAVITY_METRIC_BUFFER_NM added in Stage 1).
            ptxt = read(joinpath(leaf, "params.json"), String)
            @test occursin("\"K_LIST\"", ptxt)
            @test occursin("\"ALPHA_NM\"", ptxt)
            @test occursin("\"truncated_sides\"", ptxt)
            @test occursin("\"METHOD\"", ptxt)
            @test occursin("\"outer_polygon\"", ptxt)
            @test occursin("\"CONCAVITY_METRIC_BUFFER_NM\"", ptxt)

            # manifest.json reports the bumped loop_diagnostics_csv
            # schema_version (2 since Stage 1).
            mtxt2 = read(joinpath(leaf, "manifest.json"), String)
            @test occursin("\"loop_diagnostics_csv\"", mtxt2)
            # Schema version 2 must appear in the file (it is the only
            # schema_version: 2 entry under artifacts).
            @test occursin("\"schema_version\": 2", mtxt2)
        end
    end

    @testset "synthetic crescent: v1 misclassifies bay" begin
        # Crescent = (large disk) minus (smaller offset disk). The bay is
        # the missing-disk region; emitters do not exist there. v1's
        # alpha-shape outer polygon will SMOOTH ACROSS the concave bay
        # mouth, then point-in-polygon will say "inside" for any test
        # location in the bay. The concavity metric should run cleanly on
        # boundary-proximal interior emitters near the bay opening.
        rng = Random.MersenneTwister(7)
        cx_big, cy_big, R_big = 5.0, 5.0, 3.0
        cx_sm,  cy_sm,  R_sm  = 7.5, 5.0, 2.2
        npts = 30000
        xs = Float64[]; ys = Float64[]
        while length(xs) < npts
            xx = cx_big - R_big + 2*R_big*rand(rng)
            yy = cy_big - R_big + 2*R_big*rand(rng)
            in_big = (xx-cx_big)^2 + (yy-cy_big)^2 <= R_big^2
            in_small = (xx-cx_sm)^2 + (yy-cy_sm)^2 <= R_sm^2
            if in_big && !in_small
                push!(xs, xx); push!(ys, yy)
            end
        end

        result = classify_emitters(xs, ys; fov_um = (0.0, 10.0, 0.0, 10.0),
            params = EdgeClassifyParams(K_LIST = [8, 32], RHO_K_THRESH = 50.0,
                                        ALPHA_NM = 400.0, MEMBRANE_NM = 100.0,
                                        REFLECT_RADIUS_NM = 200.0))

        @test sum(==("outside"), result.class) +
              sum(==("membrane"), result.class) +
              sum(==("interior"), result.class) == length(xs)

        cm = compute_concavity_metric(result, xs, ys;
                                      asym_R_nm = 800.0, asym_gate = 0.10,
                                      rho_lo = 1e6)
        @test cm.n_interior >= 0
        @test cm.n_eligible >= 0
        @test cm.n_suspect >= 0
        @test length(cm.suspect_x_um) == cm.n_suspect
        @test length(cm.suspect_y_um) == cm.n_suspect
        @test length(cm.suspect_is_fov_edge) == cm.n_suspect
        @test cm.n_suspect_interior_fov + cm.n_suspect_fov_edge == cm.n_suspect
    end

    @testset "concavity metric: length validation" begin
        rng = Random.MersenneTwister(99)
        npts = 500
        u = rand(rng, npts); v = rand(rng, npts)
        r = sqrt.(u); θ = 2π .* v
        x = 5.0 .+ 2.0 .* r .* cos.(θ); y = 5.0 .+ 2.0 .* r .* sin.(θ)
        result = classify_emitters(x, y; fov_um = (0.0, 10.0, 0.0, 10.0),
            params = EdgeClassifyParams(K_LIST = [4, 16], RHO_K_THRESH = 5.0,
                                        ALPHA_NM = 800.0, REFLECT_RADIUS_NM = 100.0))
        @test_throws ArgumentError compute_concavity_metric(result, x[1:end-1], y)
        @test_throws ArgumentError compute_concavity_metric(result, x, y[1:end-1])
    end

end
