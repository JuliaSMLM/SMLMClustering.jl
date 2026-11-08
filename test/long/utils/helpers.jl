# Test helpers shared by the Core files (test/*.jl) and the Long files (test/long/*.jl).
# Lives in Long's helper folder because test/ itself holds no helper folder (lab test layout).
using SMLMData, Random

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
