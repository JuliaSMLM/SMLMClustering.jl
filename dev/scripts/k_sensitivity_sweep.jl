# kNN k-sensitivity sweep for the kNN-MRF density-regime clustering pipeline.
# Round 016, Priority 8.
#
# For each density_k k ∈ {5, 10, 15, 20, 30, 50, 80}, regenerate the synthetic
# A431-mimic at fixed 2× density ratio + patch_scale=0.75 (nominal mean rect
# length ≈ 1.5 μm — middle of the V16 plateau), run the kNN-MRF backend with
# `MRFDensityClusterConfig(n_regimes=2, density_estimator=:knn, density_k=k)`,
# record per-emitter accuracy/precision/recall, and render an SMLMRender
# CircleRender categorical TP/TN/FP/FN circle plot per k.
#
# The motivating bound is V12: kNN ball radius r_k ≈ √(k / πρ) must be smaller
# than the structure half-width or GMM regime split flips. At ρ_high=1000/μm²
# (the patch-interior density), r_k ≈ {40, 56, 69, 80, 98, 126, 160} nm for
# k ∈ {5, 10, 15, 20, 30, 50, 80}. patch_scale=0.75 gives rect half-width
# ∈ [37.5, 75] nm and ellipse semi-minor ∈ [187, 750] nm; the bound predicts
# a soft floor near k where r_k ≈ structure half-width — for the rect-thin tail
# the bound bites at k ≈ 5-10, for the ellipse bulk it stays comfortable across
# the whole sweep. Headline expectation: peak at moderate k (15-30), softening
# at very low k (sub-sampling noise → σ_log = 1/√k blows up) and at very high
# k (kNN ball spills out of patches into background → density estimate
# regresses to the global mean).
#
# Outputs in dev/scripts/output/:
#   k_sensitivity_sweep.csv                    per-k metrics
#   k_sensitivity_sweep.png                    accuracy curve, log-x in k
#   round_016_smlmrender_categorical.png       7-panel SMLMRender grid (canonical)

using Pkg
Pkg.activate(@__DIR__)

using SMLMClustering
using SMLMData
using SMLMRender
using CairoMakie
using Statistics
using Random
using Printf

# CircleRender accumulates intensity where emitter circles overlap, so individual
# RGB channels can exceed 1.0. CairoMakie's PNG backend encodes via N0f8 (8-bit
# fixed-point [0,1]) and rejects out-of-range floats. Clamp each channel by
# rebuilding the same RGB type via reflective construction (no ColorTypes dep).
clamp_rgb(c) = typeof(c)(clamp(c.r, 0.0, 1.0), clamp(c.g, 0.0, 1.0), clamp(c.b, 0.0, 1.0))

include(joinpath(@__DIR__, "simulate_a431_mimic.jl"))

const OUT_DIR = joinpath(@__DIR__, "output")

# Sweep range — log-spaced enough to see the structure, dense enough around
# the V12-default k=20 to pick up the local optimum cleanly.
const KS = [5, 10, 15, 20, 30, 50, 80]

# Patch scale fixes structure size at nominal mean rect length ≈ 1.5 μm.
# Density ratio fixed at 2× (the V13 / V16 reference contrast).
const PATCH_SCALE = 0.75
const NOMINAL_PATCH_UM = 2 * PATCH_SCALE   # ≈ 1.5 μm

# ---------------------------------------------------------------------------
# Backend predictor — Vector{Int} (1=low, 2=high) aligned to emitter order.
# ---------------------------------------------------------------------------

function predict_mrf_knn(smld, k::Int)
    cfg = MRFDensityClusterConfig(n_regimes = 2, density_estimator = :knn, density_k = k)
    (smld_out, _) = cluster(smld, cfg)
    return smld_out.metadata["mrf_regime_per_emitter"]::Vector{Int}
end

# ---------------------------------------------------------------------------
# Metrics — TP/TN/FP/FN, accuracy, precision, recall (high=positive class).
# ---------------------------------------------------------------------------
function compute_metrics(gt::Vector{Int}, pred::Vector{Int})
    n = length(gt)
    tp = tn = fp = fn = 0
    for i in 1:n
        if gt[i] == 2 && pred[i] == 2
            tp += 1
        elseif gt[i] == 1 && pred[i] == 1
            tn += 1
        elseif gt[i] == 1 && pred[i] == 2
            fp += 1
        else
            fn += 1
        end
    end
    acc = (tp + tn) / n
    precision = tp == 0 ? 0.0 : tp / (tp + fp)
    recall    = tp == 0 ? 0.0 : tp / (tp + fn)
    return (acc = acc, precision = precision, recall = recall,
            tp = tp, tn = tn, fp = fp, fn = fn)
end

# ---------------------------------------------------------------------------
# Generate the synthetic ONCE — patch_scale + ratio fixed; only k varies.
# ---------------------------------------------------------------------------
println("[sweep] simulating synthetic A431-mimic — ratio=2.0×  patch_scale=$PATCH_SCALE  (≈ $NOMINAL_PATCH_UM μm patches)")
sim_out = simulate_dataset(patch_scale = PATCH_SCALE, verbose = false)
smld    = sim_out.smld
gt      = sim_out.ground_truth
println("  n_total=$(sim_out.stats.n_total) n_low=$(sim_out.stats.n_low) n_high=$(sim_out.stats.n_high) " *
        "high_frac=$(round(100 * sim_out.stats.n_high / sim_out.stats.n_total, digits=1))%  " *
        "n_patches=$(sim_out.stats.n_patches)")

# ---------------------------------------------------------------------------
# Sweep — same SMLD, vary k.
# ---------------------------------------------------------------------------
results = NamedTuple[]
preds_by_k = Dict{Int, Vector{Int}}()

for k in KS
    println("\n[sweep] k=$k")
    t0 = time()
    pred = predict_mrf_knn(smld, k)
    elapsed_ms = round((time() - t0) * 1000, digits = 1)
    m = compute_metrics(gt, pred)
    # kNN ball radius at ρ_high = 1000 emit/μm² = 0.001 emit/nm² → r_k in nm.
    # ρ_high = rho_low + rho_high_bonus = 500 + 500 = 1000 emit/μm² = 1e-3 nm⁻²
    rho_high_nm2 = 1e-3
    r_k_nm = sqrt(k / (π * rho_high_nm2))
    push!(results, (k = k, m..., r_k_nm = r_k_nm))
    preds_by_k[k] = pred
    @printf("  k=%-3d  acc=%6.2f%% prec=%6.2f%% rec=%6.2f%%   TP=%d FP=%d FN=%d   r_k≈%.0f nm   %.0f ms\n",
            k, 100*m.acc, 100*m.precision, 100*m.recall, m.tp, m.fp, m.fn, r_k_nm, elapsed_ms)
end

# ---------------------------------------------------------------------------
# CSV
# ---------------------------------------------------------------------------
csv_path = joinpath(OUT_DIR, "k_sensitivity_sweep.csv")
open(csv_path, "w") do io
    println(io, "k,r_k_nm,n_total,accuracy,precision,recall,tp,tn,fp,fn")
    for r in sort(results, by = x -> x.k)
        n = r.tp + r.tn + r.fp + r.fn
        @printf(io, "%d,%.1f,%d,%.4f,%.4f,%.4f,%d,%d,%d,%d\n",
                r.k, r.r_k_nm, n, r.acc, r.precision, r.recall,
                r.tp, r.tn, r.fp, r.fn)
    end
end
println("\n[sweep] wrote $csv_path")

# ---------------------------------------------------------------------------
# Accuracy curve — kNN-MRF across the k sweep, log-x.
# ---------------------------------------------------------------------------
fig = Figure(size = (900, 520))
ax  = Axis(fig[1, 1];
           title  = @sprintf("kNN-MRF k-sensitivity sweep (5×5 μm A431-mimic, 2× ratio, patch_scale=%.2f ≈ %.1f μm patches)",
                              PATCH_SCALE, NOMINAL_PATCH_UM),
           xlabel = "density_k (log scale)",
           ylabel = "accuracy",
           xscale = log10,
           xticks = (KS, [string(k) for k in KS]))
hlines!(ax, [0.75]; color = (:gray60, 0.7), linestyle = :dash)
hlines!(ax, [0.85]; color = (:gray40, 0.7), linestyle = :dot)
text!(ax, KS[1] * 0.9, 0.755; text = "75% gate", color = :gray40, fontsize = 11)
text!(ax, KS[1] * 0.9, 0.855; text = "85% gate", color = :gray30, fontsize = 11)
ylims!(ax, 0.4, 1.0)

rows = sort(results, by = x -> x.k)
xs_pts = [r.k   for r in rows]
accs   = [r.acc for r in rows]
lines!(ax,    xs_pts, accs; color = :crimson, linewidth = 2.5, label = "mrf_knn")
scatter!(ax,  xs_pts, accs; color = :crimson, markersize = 12, strokewidth = 0)
# annotate r_k at each point
for r in rows
    text!(ax, r.k, r.acc + 0.025;
          text = @sprintf("r=%.0f nm", r.r_k_nm),
          align = (:center, :bottom), fontsize = 10, color = (:crimson, 0.7))
end
axislegend(ax; position = :rb, framevisible = false)
save(joinpath(OUT_DIR, "k_sensitivity_sweep.png"), fig; px_per_unit = 2)
println("[sweep] wrote k_sensitivity_sweep.png")

# ---------------------------------------------------------------------------
# SMLMRender categorical-color circle plot — canonical Round 016 deliverable.
# 7 panels, one per k value, kNN-MRF TP/TN/FP/FN per emitter.
# ---------------------------------------------------------------------------

# Encode TP/TN/FP/FN as integer ids 1..4 for SMLMRender categorical mapping.
#   1 = TN    2 = TP    3 = FP    4 = FN
function category_per_emitter(gt::Vector{Int}, pred::Vector{Int})
    cats = Vector{Int}(undef, length(gt))
    for i in eachindex(gt)
        if gt[i] == 1 && pred[i] == 1
            cats[i] = 1
        elseif gt[i] == 2 && pred[i] == 2
            cats[i] = 2
        elseif gt[i] == 1 && pred[i] == 2
            cats[i] = 3
        else
            cats[i] = 4
        end
    end
    return cats
end

function render_per_k_panel(smld, gt, pred; zoom::Int = 8)
    cats = category_per_emitter(gt, pred)
    smld_cat = deepcopy(smld)
    for (i, em) in enumerate(smld_cat.emitters)
        em.id = cats[i]
    end
    (img, _) = render(smld_cat;
                      strategy = CircleRender(),
                      color_by = :id,
                      categorical = true,
                      colormap = :tab10,
                      zoom = zoom)
    return clamp_rgb.(img)
end

println("\n[smlmrender] building categorical panel grid…")
panel_imgs = Array{Any}(undef, length(KS))
for (idx, k) in enumerate(KS)
    pred = preds_by_k[k]
    panel_imgs[idx] = render_per_k_panel(smld, gt, pred; zoom = 8)
    println("  k=$k  rendered  size=$(size(panel_imgs[idx]))")
end

# Compose into a 1xN multipanel figure via CairoMakie image! per axis.
fig2 = Figure(size = (260 * length(KS) + 60, 460))
Label(fig2[0, 1:length(KS)],
      @sprintf("kNN-MRF per-emitter classification — k-sensitivity sweep (2× ratio, ~%.1f μm patches; SMLMRender CircleRender, categorical id 1=TN 2=TP 3=FP 4=FN)",
                NOMINAL_PATCH_UM);
      fontsize = 14)
for (idx, k) in enumerate(KS)
    pred = preds_by_k[k]
    m = compute_metrics(gt, pred)
    img = panel_imgs[idx]
    ax = Axis(fig2[1, idx];
              title = @sprintf("k=%d   acc=%.1f%%   FN=%d  FP=%d",
                                k, 100*m.acc, m.fn, m.fp),
              aspect = DataAspect())
    hidedecorations!(ax)
    image!(ax, rotr90(img))
end
save(joinpath(OUT_DIR, "round_016_smlmrender_categorical.png"), fig2; px_per_unit = 2)
println("[smlmrender] wrote round_016_smlmrender_categorical.png")

println("\n[sweep] DONE — outputs in $OUT_DIR")
for f in sort(filter(x -> startswith(x, "k_sensitivity_sweep") || startswith(x, "round_016"),
                     readdir(OUT_DIR)))
    sz = stat(joinpath(OUT_DIR, f)).size
    println("  $f  ($(round(sz/1024, digits=1)) KiB)")
end
