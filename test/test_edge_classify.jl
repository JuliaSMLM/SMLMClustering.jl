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

            # loop_diagnostics.csv columns + heuristic_type LAST
            lcsv = readlines(joinpath(leaf, "loop_diagnostics.csv"))
            cols = split(lcsv[1], ",")
            @test cols == ["loop_id","vertex_count","area_um2","n_emitters_inside",
                           "frac_in_fov","frac_dense","median_rhoK","heuristic_type"]
            @test cols[end] == "heuristic_type"

            # manifest.json basic well-formedness
            mtxt = read(joinpath(leaf, "manifest.json"), String)
            @test occursin("\"schema_version\"", mtxt)
            @test occursin("\"classified_tsv\"", mtxt)
            @test occursin("\"manifest.json\"", mtxt) == false   # not self-listed
            @test occursin("\"out_dir\"", mtxt)
            @test occursin("\"leaf_dir\"", mtxt)

            # params.json carries the actual params used
            ptxt = read(joinpath(leaf, "params.json"), String)
            @test occursin("\"K_LIST\"", ptxt)
            @test occursin("\"ALPHA_NM\"", ptxt)
            @test occursin("\"truncated_sides\"", ptxt)
        end
    end

end
