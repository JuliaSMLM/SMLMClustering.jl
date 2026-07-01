using Test
using SMLMClustering

@testset "Region (MultiCellMask)" begin
    sq(x0, y0, s) = [(Float64(x0), Float64(y0)), (Float64(x0 + s), Float64(y0)),
                     (Float64(x0 + s), Float64(y0 + s)), (Float64(x0), Float64(y0 + s))]

    @testset "simple cell" begin
        m = build_mask([sq(0, 0, 10)])
        @test length(m) == 1
        @test m[1] isa CellPolygon
        @test isempty(m[1].holes)
        @test in_region(5.0, 5.0, m)
        @test !in_region(20.0, 20.0, m)
        @test isapprox(region_area(m), 100.0; atol = 1e-6)
    end

    @testset "hole keep vs discard" begin
        outer = sq(0, 0, 10); hole = sq(3, 3, 4)
        mk = build_mask([outer, hole]; keep_internal = true)
        @test length(mk) == 1
        @test length(mk[1].holes) == 1
        @test !in_region(5.0, 5.0, mk)              # (5,5) is in the hole
        @test in_region(1.0, 1.0, mk)               # inside cell, not in hole
        @test isapprox(region_area(mk), 100.0 - 16.0; atol = 1e-6)
        md = build_mask([outer, hole]; keep_internal = false)
        @test isempty(md[1].holes)
        @test in_region(5.0, 5.0, md)               # hole discarded → solid
        @test isapprox(region_area(md), 100.0; atol = 1e-6)
    end

    @testset "multiple cells + debris cutoff" begin
        big = sq(0, 0, 10); other = sq(20, 0, 10); tiny = sq(40, 40, 1)
        @test length(build_mask([big, other])) == 2                  # equal area → both kept
        @test length(build_mask([big, tiny])) == 1                   # tiny < 1/3 → dropped
        @test length(build_mask([big, tiny]; min_cell_frac = 0)) == 2  # 0 keeps everything
    end

    @testset "self-touching loop split" begin
        # figure-8: the shared vertex (5,5) splits the loop into two simple rings
        fig8 = [(0.0, 0.0), (10.0, 0.0), (5.0, 5.0), (10.0, 10.0), (0.0, 10.0), (5.0, 5.0)]
        m = build_mask([fig8]; min_cell_frac = 0)
        @test length(m) >= 1
        @test all(c -> length(c.outer) >= 3, m)
    end

    @testset "empty input" begin
        @test isempty(build_mask(Vector{NTuple{2,Float64}}[]))
        @test !in_region(0.0, 0.0, CellPolygon[])
        @test region_area(CellPolygon[]) == 0.0
    end
end
