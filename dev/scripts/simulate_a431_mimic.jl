# Build a synthetic A431-mimic SMLM dataset.
#
# Geometry:
#   5x5 um total ROI; rightmost 1/8 (x in [4.375, 5.0]) is empty (cell edge).
#   Active area = [0, 4.375] x [0, 5.0] = 21.875 um^2.
#   25% of active area in high-density patches:
#     - 60% elongated rectangles, length 1-3 um, aspect ratio 5-20
#     - 40% rotated ellipses, semi-major 0.5-1.0 um, semi-minor in [0.5, semi_major]
#   Densities:
#     low  = 500 emit/um^2
#     high = 1000 emit/um^2  (background 500 + patch contribution 500 = 1000 inside patches)
#
# Outputs (saved as JLD2 in dev/scripts/output/):
#   synthetic_smld.jld2   — BasicSMLD + sidecar Vector{Int} ground-truth labels (1=low, 2=high)
#                          + Vector{Patch} patch geometry list (for renderer overlays)

using Pkg
Pkg.activate(@__DIR__)

using Random
using SMLMData
using JLD2

const ROI_SIZE_UM   = 5.0
const EDGE_CUT_FRAC = 1/8
const ACTIVE_X_MAX  = ROI_SIZE_UM * (1 - EDGE_CUT_FRAC)  # 4.375
const ACTIVE_Y_MAX  = ROI_SIZE_UM
const ACTIVE_AREA   = ACTIVE_X_MAX * ACTIVE_Y_MAX        # 21.875 um^2
const TARGET_HIGH_FRAC = 0.25
const RHO_LOW   = 500.0   # emitters per um^2
const RHO_HIGH_BONUS = 500.0  # added on top of low inside patches → total 1000 in patches
const SIM_SEED  = 20260429

# ---------------------------------------------------------------------------
# Patch geometry (rectangle or ellipse), with bbox + point-in-shape predicate.
# ---------------------------------------------------------------------------
struct Patch
    kind::Symbol      # :rect or :ellipse
    cx::Float64
    cy::Float64
    a::Float64        # rect: half-length along major axis; ellipse: semi-major
    b::Float64        # rect: half-width;                   ellipse: semi-minor
    theta::Float64    # rotation (radians)
    area::Float64
end

function patch_bbox(p::Patch)
    # Conservative axis-aligned bounding box of the rotated shape.
    c = abs(cos(p.theta)); s = abs(sin(p.theta))
    dx = p.a * c + p.b * s
    dy = p.a * s + p.b * c
    return (p.cx - dx, p.cy - dy, p.cx + dx, p.cy + dy)
end

function in_patch(p::Patch, x::Float64, y::Float64)
    # Transform to patch-local frame.
    dx = x - p.cx
    dy = y - p.cy
    c = cos(-p.theta); s = sin(-p.theta)
    xl = c * dx - s * dy
    yl = s * dx + c * dy
    if p.kind === :rect
        return abs(xl) <= p.a && abs(yl) <= p.b
    else  # :ellipse
        return (xl / p.a)^2 + (yl / p.b)^2 <= 1.0
    end
end

# Reject placement if patch's bbox extends outside the active area, or any
# already-placed patch shape (NOT just its bbox) overlaps with the candidate.
# Shape-shape overlap is approximated by sampling 200 points uniformly in the
# candidate's bbox and checking if any is in BOTH the candidate AND any
# existing patch. Faster than SAT for mixed rect/ellipse shapes and uniformly
# accurate enough for placement at our scale.
function placement_ok(p::Patch, placed::Vector{Patch}, rng)
    (xlo, ylo, xhi, yhi) = patch_bbox(p)
    (xlo >= 0 && ylo >= 0 && xhi <= ACTIVE_X_MAX && yhi <= ACTIVE_Y_MAX) || return false
    isempty(placed) && return true
    # Cheap reject: if bboxes don't overlap any existing, accept.
    for q in placed
        (qxlo, qylo, qxhi, qyhi) = patch_bbox(q)
        bbox_overlap = !(xhi < qxlo || xlo > qxhi || yhi < qylo || ylo > qyhi)
        bbox_overlap || continue
        # bbox overlaps q → check shape overlap by sampling.
        for _ in 1:200
            x = xlo + rand(rng) * (xhi - xlo)
            y = ylo + rand(rng) * (yhi - ylo)
            if in_patch(p, x, y) && in_patch(q, x, y)
                return false  # shape overlap detected
            end
        end
        # 200 samples found no overlap → treat as non-overlapping
    end
    return true
end

# Draw one patch (kind sampled), random angle, random center inside active area.
function sample_patch(rng, kind::Symbol)
    if kind === :rect
        L = rand(rng) * 2.0 + 1.0           # length in [1, 3]
        ar = rand(rng) * 15.0 + 5.0         # aspect ratio in [5, 20]
        W = L / ar
        a = L / 2; b = W / 2
        area = L * W
    else
        a = rand(rng) * 0.5 + 0.5           # semi-major in [0.5, 1.0]
        b = a * (rand(rng) * 0.5 + 0.5)     # semi-minor in [0.5*a, a]
        area = π * a * b
    end
    theta = rand(rng) * 2π
    # cx, cy will be re-rolled in the placement loop.
    return Patch(kind, 0.0, 0.0, a, b, theta, area)
end

function place_patches(rng)
    placed = Patch[]
    placed_area = 0.0
    target_area = TARGET_HIGH_FRAC * ACTIVE_AREA
    failed_attempts = 0
    while placed_area < target_area
        kind = rand(rng) < 0.6 ? :rect : :ellipse
        p_template = sample_patch(rng, kind)
        success = false
        for _ in 1:100
            # Roll center until bbox fits in active area.
            cx = rand(rng) * ACTIVE_X_MAX
            cy = rand(rng) * ACTIVE_Y_MAX
            p = Patch(p_template.kind, cx, cy, p_template.a, p_template.b,
                     p_template.theta, p_template.area)
            if placement_ok(p, placed, rng)
                push!(placed, p)
                placed_area += p.area
                success = true
                break
            end
        end
        if !success
            failed_attempts += 1
            if failed_attempts > 50
                @info "Stopped placing patches after 50 consecutive failures" placed_area target_area
                break
            end
        else
            failed_attempts = 0
        end
    end
    return placed, placed_area
end

# ---------------------------------------------------------------------------
# Emitter generation. Background (low) + per-patch additions (bringing patches
# up to high density). Returns parallel vectors of (x, y, ground_truth_label).
# ---------------------------------------------------------------------------
function poisson_count(rng, mean_count::Float64)
    # Knuth for small λ; Gaussian approximation N(λ, λ) for λ ≥ 30 (Knuth
    # underflows exp(-λ) to 0 at large λ and never terminates correctly).
    if mean_count < 30
        L = exp(-mean_count)
        k = 0
        p = 1.0
        while true
            k += 1
            p *= rand(rng)
            p <= L && return k - 1
        end
    else
        return max(0, round(Int, mean_count + sqrt(mean_count) * randn(rng)))
    end
end

function generate_emitters(rng, patches; rho_low::Float64 = RHO_LOW,
                            rho_high_bonus::Float64 = RHO_HIGH_BONUS)
    xs = Float64[]; ys = Float64[]; gt = Int[]

    # 1. Background: Poisson over the full active area at rho_low.
    n_bg = poisson_count(rng, rho_low * ACTIVE_AREA)
    for _ in 1:n_bg
        x = rand(rng) * ACTIVE_X_MAX
        y = rand(rng) * ACTIVE_Y_MAX
        push!(xs, x); push!(ys, y)
        # Label by whether this point happens to land in a patch.
        in_high = false
        for p in patches
            if in_patch(p, x, y); in_high = true; break; end
        end
        push!(gt, in_high ? 2 : 1)
    end

    # 2. Per-patch additions: Poisson(rho_high_bonus * patch.area), rejection-sampled
    # uniformly inside the patch shape using the bbox.
    for p in patches
        n_add = poisson_count(rng, rho_high_bonus * p.area)
        (xlo, ylo, xhi, yhi) = patch_bbox(p)
        added = 0
        while added < n_add
            x = xlo + rand(rng) * (xhi - xlo)
            y = ylo + rand(rng) * (yhi - ylo)
            if in_patch(p, x, y)
                push!(xs, x); push!(ys, y); push!(gt, 2)
                added += 1
            end
        end
    end

    return xs, ys, gt
end

# ---------------------------------------------------------------------------
# K distribution per @genmab handwave: 80% K=1; 15% from Geometric(0.5)+1
# truncated to [2,5]; 5% from Geometric(0.1)+5 truncated to [6,50].
# Photons = 1000 * K.
# ---------------------------------------------------------------------------
function sample_k(rng)
    r = rand(rng)
    if r < 0.80
        return 1
    elseif r < 0.95
        for _ in 1:50
            k = 1 + ceil(Int, log(rand(rng)) / log(0.5))   # Geometric(0.5)+1
            2 <= k <= 5 && return k
        end
        return 3
    else
        for _ in 1:50
            k = 5 + ceil(Int, log(rand(rng)) / log(0.9))   # Geometric(0.1)+5
            6 <= k <= 50 && return k
        end
        return 10
    end
end

function build_smld(xs::Vector{Float64}, ys::Vector{Float64}, rng;
                    rho_low::Float64 = RHO_LOW,
                    rho_high_bonus::Float64 = RHO_HIGH_BONUS,
                    seed::Int = SIM_SEED)
    cam = SMLMData.IdealCamera(1:64, 1:64, 0.1)
    n = length(xs)
    emitters = Vector{SMLMData.Emitter2DFit{Float64}}(undef, n)
    for i in 1:n
        K = sample_k(rng)
        photons   = 1000.0 * K
        emitters[i] = SMLMData.Emitter2DFit{Float64}(
            xs[i], ys[i],          # x, y (μm)
            photons, 10.0,         # photons, bg
            0.008, 0.008,          # σ_x, σ_y (μm = 8 nm)
            50.0, 0.5;             # σ_photons, σ_bg
            frame = 1, dataset = 1,
        )
    end
    return SMLMData.BasicSMLD(emitters, cam, 1, 1,
                              Dict{String,Any}(
                                  "simulation_seed" => seed,
                                  "roi_size_um" => ROI_SIZE_UM,
                                  "edge_cut_frac" => EDGE_CUT_FRAC,
                                  "rho_low" => rho_low,
                                  "rho_high" => rho_low + rho_high_bonus,
                                  "target_high_area_frac" => TARGET_HIGH_FRAC,
                              ))
end

# ---------------------------------------------------------------------------
# Reusable entry point — same patch geometry across calls when seed fixed.
# Returns NamedTuple for consumers (sweep scripts, calibration scripts).
# ---------------------------------------------------------------------------
function simulate_dataset(; rho_low::Float64 = RHO_LOW,
                          rho_high_bonus::Float64 = RHO_HIGH_BONUS,
                          seed::Int = SIM_SEED,
                          verbose::Bool = false)
    rng = Xoshiro(seed)
    verbose && println("[simulate] placing patches…")
    patches, placed_area = place_patches(rng)
    actual_high_frac = placed_area / ACTIVE_AREA
    verbose && println("  placed $(length(patches)) patches; cumulative area = $(round(placed_area, digits=3)) μm² ",
                       "($(round(100*actual_high_frac, digits=1))% of active area)")
    verbose && println("[simulate] generating emitters (rho_low=$rho_low, rho_high_bonus=$rho_high_bonus)…")
    xs, ys, gt = generate_emitters(rng, patches; rho_low = rho_low, rho_high_bonus = rho_high_bonus)
    n_total = length(xs); n_low = count(==(1), gt); n_high = count(==(2), gt)
    verbose && println("  total = $n_total | low = $n_low | high = $n_high | high frac = $(round(100*n_high/n_total, digits=1))%")
    smld = build_smld(xs, ys, rng; rho_low = rho_low, rho_high_bonus = rho_high_bonus, seed = seed)
    return (smld = smld, ground_truth = gt, patches = patches,
            stats = (n_total = n_total, n_low = n_low, n_high = n_high,
                     n_patches = length(patches),
                     actual_high_area_frac = actual_high_frac,
                     rho_low = rho_low, rho_high_bonus = rho_high_bonus,
                     density_ratio = (rho_low + rho_high_bonus) / rho_low))
end

# ---------------------------------------------------------------------------
# Main — only runs when this file is executed as a script.
# Other scripts (sweeps, calibration demos) `include()` this file to get
# `simulate_dataset(...)` and the patch-geometry helpers without re-running
# the default-density simulate-and-save here.
# ---------------------------------------------------------------------------
function main()
    out = simulate_dataset(verbose = true)
    out_dir = joinpath(@__DIR__, "output")
    isdir(out_dir) || mkpath(out_dir)
    out_path = joinpath(out_dir, "synthetic_smld.jld2")
    jldsave(out_path;
            smld = out.smld,
            ground_truth = out.ground_truth,
            patches = out.patches,
            stats = out.stats)
    println("[simulate] saved $out_path")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
