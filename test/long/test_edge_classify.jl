# EdgeClassify on SMLMSim-simulated cells: a dense square cell in a sparse background,
# with blinking and localization noise from SMLMSim. The structural checks of the
# thorough tier, plus ground truth that only a simulation provides.
using Test
using Random
using SMLMClustering
using SMLMData
using SMLMSim

# One emitter per site, placed uniformly by SMLMSim at density ρ_cell (µm⁻²) inside a
# square cell of side `cell` centered in a `field` × `field` µm FOV, and at ρ_bg over the
# whole FOV. SMLMSim then blinks every emitter over `nframes` frames and adds
# localization noise. Each localization's track_id is its true emitter's index.
function _sim_cell(seed; ρ_cell, ρ_bg, cell = 6.0, field = 10.0, nframes = 2000)
    Random.seed!(seed)
    site = Nmer2D(n = 1, d = 0.0)
    off = (field - cell) / 2
    cx, cy, _ = uniform2D(ρ_cell, site, cell, cell)
    bx, by, _ = uniform2D(ρ_bg, site, field, field)
    x = vcat(cx .+ off, bx); y = vcat(cy .+ off, by)
    em = [Emitter2DFit{Float64}(x[i], y[i], 1.0, 0.0, 0.0, 0.0, 0.0, 0.0;
                                frame = 1, dataset = 1, track_id = i, id = i) for i in eachindex(x)]
    npx = round(Int, field / 0.1)
    smld, _ = simulate(StaticSMLMConfig(ndatasets = 1, nframes = nframes);
                       starting_conditions = em, camera = IdealCamera(1:npx, 1:npx, 0.1))
    return smld, length(cx), off, off + cell
end

@testset "EdgeClassify on SMLMSim-simulated cells" begin

    @testset "dense cell at $(Int(ρc)) emitters/µm²" for (seed, ρc, ρb) in
            ((1, 400.0, 4.0), (2, 200.0, 4.0), (3, 800.0, 8.0))
        smld, ncell, lo, hi = _sim_cell(seed; ρ_cell = ρc, ρ_bg = ρb)
        _, info = classify_emitters(smld, KdeValleyConfig())

        # Structural checks.
        @test info.n_outside + info.n_membrane + info.n_interior == info.n_emitters
        @test info.n_emitters == length(smld.emitters)
        @test length(info.cells) == 1
        @test info.outer_polygon == info.cells[1].outer
        @test all((info.class .!= :outside) .== info.inside_outer)
        @test 0.8 < interior_fraction(info) < 0.995

        # Ground truth: localizations of cell emitters are labeled in the cell; background
        # localizations at least 0.5 µm outside the square are labeled outside.
        incell = [e.track_id <= ncell for e in smld.emitters]
        far = [!incell[k] && (min(e.x, e.y) < lo - 0.5 || max(e.x, e.y) > hi + 0.5)
               for (k, e) in enumerate(smld.emitters)]
        @test count(incell .& (info.class .!= :outside)) / count(incell) > 0.95
        @test count(far .& (info.class .== :outside)) / count(far) > 0.98
    end

    @testset "background only: completes, counts add up" begin
        Random.seed!(9)
        bx, by, _ = uniform2D(4.0, Nmer2D(n = 1, d = 0.0), 10.0, 10.0)
        em = [Emitter2DFit{Float64}(bx[i], by[i], 1.0, 0.0, 0.0, 0.0, 0.0, 0.0;
                                    frame = 1, dataset = 1, track_id = i, id = i) for i in eachindex(bx)]
        smld, _ = simulate(StaticSMLMConfig(ndatasets = 1, nframes = 2000);
                           starting_conditions = em, camera = IdealCamera(1:100, 1:100, 0.1))
        _, info = classify_emitters(smld, KdeValleyConfig())
        @test info.n_outside + info.n_membrane + info.n_interior == info.n_emitters
        @test info.n_emitters == length(smld.emitters)
    end

end
