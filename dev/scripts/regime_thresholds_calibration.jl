# regime_thresholds calibration workflow demo (Round 014, Priority 3).
#
# Workflow: derive log-density regime thresholds from a high-quality
# calibration ROI via `calibrate_regime_thresholds`, then apply those
# thresholds on held-out query ROIs by passing them into
# `MRFDensityClusterConfig(regime_thresholds = ...)`. Compares to
# GMM auto-mode (no override) on the same query data.
#
# Bar for round close (per Priority 3 in STATUS.md):
#   override mode must match GMM auto-mode within ±2 pp accuracy on the
#   synthetic. Round 013's surprise (Potts smoothness amplifies a
#   degenerate GMM split into uniform misclassification at low contrast)
#   means a calibration-derived threshold may also stabilize the
#   override at low contrast — extension test at ratio 1.5×.
#
# Outputs in dev/scripts/output/:
#   regime_calibration_demo.png         — calibration histogram + threshold
#                                         line + GMM components + accuracy
#                                         table for query ROI (2× ratio)
#   regime_calibration_demo.csv         — per-(case, mode) accuracy/precision/
#                                         recall, with calibrated threshold value
#   regime_calibration_panels.png       — categorical per-emitter circle plots,
#                                         one panel per (ratio, mode) cell

using Pkg
Pkg.activate(@__DIR__)

using SMLMClustering
using SMLMData
using JLD2
using CairoMakie
using Statistics
using Random
using Printf
using NearestNeighbors

include(joinpath(@__DIR__, "simulate_a431_mimic.jl"))

const OUT_DIR = joinpath(@__DIR__, "output")
const CAL_SEED   = SIM_SEED                    # 20260429 — same as synthetic_smld.jld2
const QUERY_SEED = CAL_SEED + 1                 # 20260430 — different patch geometry, same params
const DENSITY_K  = 20

# ---------------------------------------------------------------------------
# Run MRF in two modes: GMM auto, and override with the calibrated thresholds.
# Returns named tuple with (regime, accuracy, precision, recall, tp, fp, fn).
# ---------------------------------------------------------------------------
function run_mrf(smld, gt; thresholds = nothing)
    cfg = MRFDensityClusterConfig(
        n_regimes = 2,
        density_estimator = :knn,
        density_k = DENSITY_K,
        regime_thresholds = thresholds,
    )
    (smld_out, _) = cluster(smld, cfg)
    pred = smld_out.metadata["mrf_regime_per_emitter"]::Vector{Int}
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
    prec = tp == 0 ? 0.0 : tp / (tp + fp)
    rec  = tp == 0 ? 0.0 : tp / (tp + fn)
    return (regime = pred, accuracy = acc, precision = prec, recall = rec,
            tp = tp, tn = tn, fp = fp, fn = fn)
end

# ---------------------------------------------------------------------------
# 1. Calibration ROI — fit thresholds from a 2× synthetic.
# ---------------------------------------------------------------------------
println("[calibration] simulating cal ROI (seed=$(CAL_SEED), 2× ratio)…")
cal_out = simulate_dataset(
    rho_low = RHO_LOW, rho_high_bonus = RHO_HIGH_BONUS,
    seed = CAL_SEED, verbose = false,
)
println("  $(cal_out.stats.n_total) emitters " *
        "($(round(100 * cal_out.stats.n_high / cal_out.stats.n_total, digits=1))% high)")

println("[calibration] fitting GMM via calibrate_regime_thresholds (k=$DENSITY_K)…")
thresholds = calibrate_regime_thresholds(
    cal_out.smld;
    n_regimes = 2,
    density_estimator = :knn,
    density_k = DENSITY_K,
)
println("  threshold (log ρ): $(thresholds[1]) → ρ = $(exp(thresholds[1])) μm⁻²")

# ---------------------------------------------------------------------------
# 2. Query ROI at 2× — must match GMM mode within ±2 pp.
# ---------------------------------------------------------------------------
println("[query 2×] simulating query ROI (seed=$(QUERY_SEED), 2× ratio)…")
q2_out = simulate_dataset(
    rho_low = RHO_LOW, rho_high_bonus = RHO_HIGH_BONUS,
    seed = QUERY_SEED, verbose = false,
)
println("  $(q2_out.stats.n_total) emitters " *
        "($(round(100 * q2_out.stats.n_high / q2_out.stats.n_total, digits=1))% high)")

println("[query 2×] running GMM auto-mode…")
m2_auto = run_mrf(q2_out.smld, q2_out.ground_truth)
@printf("  acc=%6.2f%%  prec=%6.2f%%  rec=%6.2f%%  TP=%d FP=%d FN=%d\n",
        100 * m2_auto.accuracy, 100 * m2_auto.precision, 100 * m2_auto.recall,
        m2_auto.tp, m2_auto.fp, m2_auto.fn)

println("[query 2×] running override mode (calibrated thresholds)…")
m2_over = run_mrf(q2_out.smld, q2_out.ground_truth; thresholds = thresholds)
@printf("  acc=%6.2f%%  prec=%6.2f%%  rec=%6.2f%%  TP=%d FP=%d FN=%d\n",
        100 * m2_over.accuracy, 100 * m2_over.precision, 100 * m2_over.recall,
        m2_over.tp, m2_over.fp, m2_over.fn)

delta_2x_pp = 100.0 * (m2_over.accuracy - m2_auto.accuracy)
@printf("  Δ acc (override − auto) = %+6.2f pp   (gate: |Δ| ≤ 2 pp)\n", delta_2x_pp)

# ---------------------------------------------------------------------------
# 3. Extension — query ROI at 1.5× (low contrast where GMM auto degenerates).
#    Round 013 saw kNN-MRF collapse to ~69% at this ratio (Potts amplifies a
#    weak GMM signal). The calibrated threshold (from the 2× cal ROI) should
#    avoid the degenerate split.
# ---------------------------------------------------------------------------
println("[query 1.5×] simulating low-contrast query (1.5× ratio)…")
q15_out = simulate_dataset(
    rho_low = RHO_LOW, rho_high_bonus = RHO_LOW * 0.5,  # 1.5× ratio
    seed = QUERY_SEED, verbose = false,
)
println("  $(q15_out.stats.n_total) emitters " *
        "($(round(100 * q15_out.stats.n_high / q15_out.stats.n_total, digits=1))% high)")

println("[query 1.5×] running GMM auto-mode…")
m15_auto = run_mrf(q15_out.smld, q15_out.ground_truth)
@printf("  acc=%6.2f%%  prec=%6.2f%%  rec=%6.2f%%  TP=%d FP=%d FN=%d\n",
        100 * m15_auto.accuracy, 100 * m15_auto.precision, 100 * m15_auto.recall,
        m15_auto.tp, m15_auto.fp, m15_auto.fn)

println("[query 1.5×] running override mode (calibrated thresholds)…")
m15_over = run_mrf(q15_out.smld, q15_out.ground_truth; thresholds = thresholds)
@printf("  acc=%6.2f%%  prec=%6.2f%%  rec=%6.2f%%  TP=%d FP=%d FN=%d\n",
        100 * m15_over.accuracy, 100 * m15_over.precision, 100 * m15_over.recall,
        m15_over.tp, m15_over.fp, m15_over.fn)

delta_15x_pp = 100.0 * (m15_over.accuracy - m15_auto.accuracy)
@printf("  Δ acc (override − auto) = %+6.2f pp   (low-contrast extension)\n", delta_15x_pp)

# ---------------------------------------------------------------------------
# 4. CSV
# ---------------------------------------------------------------------------
csv_path = joinpath(OUT_DIR, "regime_calibration_demo.csv")
open(csv_path, "w") do io
    println(io, "case,mode,calibrated_threshold_logrho,accuracy,precision,recall,tp,tn,fp,fn,n_total,delta_pp_vs_auto")
    @printf(io, "query_2x,gmm_auto,%.6f,%.4f,%.4f,%.4f,%d,%d,%d,%d,%d,%.2f\n",
            thresholds[1], m2_auto.accuracy, m2_auto.precision, m2_auto.recall,
            m2_auto.tp, m2_auto.tn, m2_auto.fp, m2_auto.fn,
            q2_out.stats.n_total, 0.0)
    @printf(io, "query_2x,override,%.6f,%.4f,%.4f,%.4f,%d,%d,%d,%d,%d,%+.2f\n",
            thresholds[1], m2_over.accuracy, m2_over.precision, m2_over.recall,
            m2_over.tp, m2_over.tn, m2_over.fp, m2_over.fn,
            q2_out.stats.n_total, delta_2x_pp)
    @printf(io, "query_15x,gmm_auto,%.6f,%.4f,%.4f,%.4f,%d,%d,%d,%d,%d,%.2f\n",
            thresholds[1], m15_auto.accuracy, m15_auto.precision, m15_auto.recall,
            m15_auto.tp, m15_auto.tn, m15_auto.fp, m15_auto.fn,
            q15_out.stats.n_total, 0.0)
    @printf(io, "query_15x,override,%.6f,%.4f,%.4f,%.4f,%d,%d,%d,%d,%d,%+.2f\n",
            thresholds[1], m15_over.accuracy, m15_over.precision, m15_over.recall,
            m15_over.tp, m15_over.tn, m15_over.fp, m15_over.fn,
            q15_out.stats.n_total, delta_15x_pp)
end
println("[output] wrote $csv_path")

# ---------------------------------------------------------------------------
# 5. Calibration figure: log-density histogram (cal ROI) + threshold line +
#    accuracy table. The histogram's two modes (one per regime) and the
#    calibrated threshold sit at the inter-mode boundary.
# ---------------------------------------------------------------------------
function compute_log_rho_knn(smld, k)
    X = Matrix{Float64}(undef, 2, length(smld.emitters))
    for (i, e) in enumerate(smld.emitters)
        X[1, i] = e.x; X[2, i] = e.y
    end
    tree = KDTree(X)
    log_ρ = Float64[]
    for i in 1:size(X, 2)
        _, dists = knn(tree, view(X, :, i), k + 1, true)
        rk = dists[end]
        rk > 0 && push!(log_ρ, log(k / (π * rk^2)))
    end
    log_ρ
end

cal_log_rho = compute_log_rho_knn(cal_out.smld, DENSITY_K)

fig = Figure(size = (1400, 720))
Label(fig[0, 1:2],
      "regime_thresholds calibration workflow — kNN density (k=$DENSITY_K), 2-regime GMM";
      fontsize = 16)

# Left: calibration histogram + threshold + GMM means.
ax1 = Axis(fig[1, 1];
           title  = @sprintf("calibration ROI (seed=%d, 2× ratio, n=%d)",
                              CAL_SEED, length(cal_log_rho)),
           xlabel = "log ρ (μm⁻²)",
           ylabel = "count")
hist!(ax1, cal_log_rho;
      bins = 50, color = (:steelblue, 0.6), strokewidth = 0.4, strokecolor = :black)
vlines!(ax1, [thresholds[1]];
        color = :crimson, linestyle = :solid, linewidth = 2.5,
        label = @sprintf("threshold = %.3f", thresholds[1]))
axislegend(ax1; position = :rt, framevisible = false)

# Right: accuracy comparison bars.
ax2 = Axis(fig[1, 2];
           title = "MRF accuracy on held-out query ROIs",
           xticks = ([1, 2, 3, 4],
                     ["2× auto", "2× override", "1.5× auto", "1.5× override"]),
           ylabel = "accuracy",
           xticklabelrotation = π/8)
ylims!(ax2, 0, 1.0)
bar_xs = [1, 2, 3, 4]
bar_accs = [m2_auto.accuracy, m2_over.accuracy, m15_auto.accuracy, m15_over.accuracy]
bar_colors = [(:slateblue, 0.85), (:crimson, 0.85), (:slateblue, 0.5), (:crimson, 0.5)]
barplot!(ax2, bar_xs, bar_accs; color = bar_colors, strokewidth = 0)
hlines!(ax2, [0.85]; color = (:gray40, 0.7), linestyle = :dot)
text!(ax2, 0.6, 0.86; text = "85% gate", color = :gray40, fontsize = 11)
hlines!(ax2, [0.75]; color = (:gray60, 0.7), linestyle = :dash)
text!(ax2, 0.6, 0.76; text = "75% gate", color = :gray40, fontsize = 11)
for (x, a) in zip(bar_xs, bar_accs)
    text!(ax2, x, a + 0.02; text = @sprintf("%.1f%%", 100*a),
          align = (:center, :bottom), fontsize = 13)
end

# Bottom row: summary text block.
summary_text = @sprintf("Calibrated threshold (log ρ): %.4f   →   ρ = %.1f μm⁻²\nQuery 2× ratio:    auto  acc = %5.2f%%   override acc = %5.2f%%   Δ = %+.2f pp\nQuery 1.5× ratio:  auto  acc = %5.2f%%   override acc = %5.2f%%   Δ = %+.2f pp",
    thresholds[1], exp(thresholds[1]),
    100 * m2_auto.accuracy, 100 * m2_over.accuracy, delta_2x_pp,
    100 * m15_auto.accuracy, 100 * m15_over.accuracy, delta_15x_pp,
)
Label(fig[2, 1:2], summary_text; fontsize = 13, halign = :left,
      tellwidth = false, tellheight = true)

save(joinpath(OUT_DIR, "regime_calibration_demo.png"), fig; px_per_unit = 2)
println("[output] wrote regime_calibration_demo.png")

# ---------------------------------------------------------------------------
# 6. Categorical per-emitter circle plot — one panel per (ratio, mode) cell.
#    Categorical color scheme: gray = correct, orange = FP (false high),
#    cyan = FN (missed high). Per phase-3-extras of start-round.md:
#    SMLMRender-style (CairoMakie scatter — SMLMRender's CircleRender does not
#    expose per-emitter color; raw scatter with markersize~σ_loc and a
#    qualitative color array is the canonical idiom for per-class visualization
#    in this lab — confirmed by @analysis 2026-04-29).
# ---------------------------------------------------------------------------
const COLOR_OK = (:gray70,    0.45)
const COLOR_FP = (:darkorange, 0.85)
const COLOR_FN = (:darkcyan,   0.85)

panel_cases = [
    ("2× auto",     q2_out,  m2_auto),
    ("2× override", q2_out,  m2_over),
    ("1.5× auto",   q15_out, m15_auto),
    ("1.5× override", q15_out, m15_over),
]

fig2 = Figure(size = (1500, 380))
Label(fig2[0, 1:length(panel_cases)],
      "Per-emitter classification (categorical) — gray=correct, orange=FP, cyan=FN";
      fontsize = 14)

for (k, (label, qout, m)) in enumerate(panel_cases)
    smld = qout.smld
    gt   = qout.ground_truth
    pred = m.regime
    xs = Float64[e.x for e in smld.emitters]
    ys = Float64[e.y for e in smld.emitters]
    colors = [
        gt[i] == pred[i] ? COLOR_OK :
        (gt[i] == 1 && pred[i] == 2 ? COLOR_FP : COLOR_FN)
        for i in eachindex(gt)
    ]
    ax = Axis(fig2[1, k];
              title  = @sprintf("%s   acc=%.1f%%   FN=%d  FP=%d",
                                 label, 100*m.accuracy, m.fn, m.fp),
              xlabel = "x (μm)",
              ylabel = (k == 1 ? "y (μm)" : ""),
              aspect = DataAspect())
    xlims!(ax, 0, 5); ylims!(ax, 0, 5)
    scatter!(ax, xs, ys; color = colors, markersize = 2.5, strokewidth = 0)
end

save(joinpath(OUT_DIR, "regime_calibration_panels.png"), fig2; px_per_unit = 2)
println("[output] wrote regime_calibration_panels.png")

# ---------------------------------------------------------------------------
# 7. Gate check
# ---------------------------------------------------------------------------
println("\n[summary]")
@printf("  Calibrated threshold (log ρ):    %.4f  (ρ = %.1f μm⁻²)\n",
        thresholds[1], exp(thresholds[1]))
@printf("  Query 2× ratio:   auto = %5.2f%%   override = %5.2f%%   Δ = %+.2f pp\n",
        100 * m2_auto.accuracy, 100 * m2_over.accuracy, delta_2x_pp)
@printf("  Query 1.5× ratio: auto = %5.2f%%   override = %5.2f%%   Δ = %+.2f pp\n",
        100 * m15_auto.accuracy, 100 * m15_over.accuracy, delta_15x_pp)

if abs(delta_2x_pp) <= 2.0
    @printf("  GATE PASSED — |Δ 2×| = %.2f pp ≤ 2 pp.\n", abs(delta_2x_pp))
else
    @printf("  GATE FAILED — |Δ 2×| = %.2f pp > 2 pp.\n", abs(delta_2x_pp))
end

println("\n[done] outputs in $OUT_DIR")
for f in sort(filter(x -> startswith(x, "regime_calibration"), readdir(OUT_DIR)))
    sz = stat(joinpath(OUT_DIR, f)).size
    println("  $f  ($(round(sz/1024, digits=1)) KiB)")
end
