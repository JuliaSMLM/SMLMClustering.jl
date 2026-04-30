# Density-ratio operational-floor sweep for the kNN-MRF density-regime
# clustering pipeline. Round 013, Priority 2.
#
# For each ratio R ∈ {1.2, 1.5, 2.0, 3.0, 5.0}: regenerate the synthetic
# A431-mimic at that density ratio (same patch geometry — fixed SIM_SEED —
# varying RHO_HIGH_BONUS so high = low * R), run 4 backends, record
# per-emitter accuracy/precision/recall, render a categorical TP/TN/FP/FN
# circle plot per ratio, and an aggregate curve plot.
#
# Outputs in dev/scripts/output/:
#   density_ratio_sweep.csv                          per (backend, ratio) metrics
#   density_ratio_sweep.png                          accuracy curve, all 4 backends
#   density_ratio_sweep_<backend>_panels.png         per-emitter circle plot, 5 ratios

using Pkg
Pkg.activate(@__DIR__)

using SMLMClustering
using SMLMData
using JLD2
using CairoMakie
using Statistics
using Random
using Printf

include(joinpath(@__DIR__, "simulate_a431_mimic.jl"))

const OUT_DIR = joinpath(@__DIR__, "output")
const RATIOS = [1.2, 1.5, 2.0, 3.0, 5.0]
# Use the constant from the included sim file; expose locally for clarity.
const SWEEP_RHO_LOW = RHO_LOW

# ---------------------------------------------------------------------------
# Backend predictors — Vector{Int} (1=low, 2=high) aligned to emitter order.
# kNN-MRF (k=20) is the post-Round-012 winner; sweep keeps k fixed and varies
# only the density ratio. k-sensitivity at fixed ratio is a follow-up.
# ---------------------------------------------------------------------------

function predict_mrf_knn(smld)
    cfg = MRFDensityClusterConfig(n_regimes = 2, density_estimator = :knn, density_k = 20)
    (smld_out, _) = cluster(smld, cfg)
    return smld_out.metadata["mrf_regime_per_emitter"]::Vector{Int}
end

function predict_hdbscan(smld)
    cfg = HDBSCANConfig(min_points = 5, knn_graph_k = 30)
    (smld_out, _) = cluster(smld, cfg)
    return [e.id == 0 ? 1 : 2 for e in smld_out.emitters]
end

function predict_dbscan(smld)
    cfg = DBSCANConfig(eps_nm = 100.0, min_points = 5)
    (smld_out, _) = cluster(smld, cfg)
    return [e.id == 0 ? 1 : 2 for e in smld_out.emitters]
end

# 2-component GMM EM on log ρ — unsmoothed Voronoi baseline.
function fit_2gmm(values::Vector{Float64}; max_iter = 200, tol = 1e-5)
    n = length(values)
    med = Statistics.median(values)
    μ = [Statistics.mean(values[values .< med]), Statistics.mean(values[values .>= med])]
    σ² = [Statistics.var(values), Statistics.var(values)]
    π_ = [0.5, 0.5]
    γ = zeros(n, 2)
    prev_ll = -Inf
    for _ in 1:max_iter
        for i in 1:n
            l1 = log(π_[1]) - 0.5*log(2π*σ²[1]) - (values[i]-μ[1])^2/(2σ²[1])
            l2 = log(π_[2]) - 0.5*log(2π*σ²[2]) - (values[i]-μ[2])^2/(2σ²[2])
            m = max(l1, l2)
            denom = m + log(exp(l1-m) + exp(l2-m))
            γ[i,1] = exp(l1 - denom); γ[i,2] = exp(l2 - denom)
        end
        for k in 1:2
            wk = max(sum(@view γ[:,k]), 1e-12)
            π_[k] = wk / n
            μ[k] = sum(γ[i,k] * values[i] for i in 1:n) / wk
            σ²[k] = max(sum(γ[i,k] * (values[i] - μ[k])^2 for i in 1:n) / wk, 1e-12)
        end
        ll = 0.0
        for i in 1:n
            l1 = log(π_[1]) - 0.5*log(2π*σ²[1]) - (values[i]-μ[1])^2/(2σ²[1])
            l2 = log(π_[2]) - 0.5*log(2π*σ²[2]) - (values[i]-μ[2])^2/(2σ²[2])
            m = max(l1, l2)
            ll += m + log(exp(l1-m) + exp(l2-m))
        end
        abs(ll - prev_ll) < tol && break
        prev_ll = ll
    end
    return μ, σ², π_, γ
end

function predict_voronoi_gmm(smld)
    cfg = VoronoiDensityConfig()
    (_, info) = cluster_statistics(smld, cfg)
    ρ = info.extras[:density_per_emitter]
    valid = .!isnan.(ρ)
    log_ρ = log.(ρ[valid])
    μ, _, _, γ = fit_2gmm(log_ρ)
    low_k = argmin(μ); high_k = 3 - low_k
    pred_valid = [γ[i, high_k] > γ[i, low_k] ? 2 : 1 for i in 1:length(log_ρ)]
    pred = ones(Int, length(ρ))
    pred[valid] = pred_valid
    return pred
end

# ---------------------------------------------------------------------------
# Metrics: TP/TN/FP/FN, accuracy, precision, recall (high=positive class).
# ---------------------------------------------------------------------------
function compute_metrics(gt::Vector{Int}, pred::Vector{Int})
    n = length(gt)
    @assert length(pred) == n
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
# Sweep
# ---------------------------------------------------------------------------
const BACKENDS = [
    ("mrf_knn",     predict_mrf_knn),
    ("dbscan",      predict_dbscan),
    ("hdbscan",     predict_hdbscan),
    ("voronoi_gmm", predict_voronoi_gmm),
]

results = Dict{String, Vector{NamedTuple}}()
for (name, _) in BACKENDS
    results[name] = NamedTuple[]
end

preds_by_ratio = Dict{Float64, Dict{String, Vector{Int}}}()
gt_by_ratio    = Dict{Float64, Vector{Int}}()
xs_by_ratio    = Dict{Float64, Vector{Float64}}()
ys_by_ratio    = Dict{Float64, Vector{Float64}}()

for ratio in RATIOS
    println("[sweep] ratio=$(ratio)x")
    bonus = SWEEP_RHO_LOW * (ratio - 1)
    out   = simulate_dataset(rho_low = SWEEP_RHO_LOW, rho_high_bonus = bonus, verbose = false)
    smld  = out.smld; gt = out.ground_truth
    xs    = Float64[e.x for e in smld.emitters]
    ys    = Float64[e.y for e in smld.emitters]

    gt_by_ratio[ratio] = gt
    xs_by_ratio[ratio] = xs
    ys_by_ratio[ratio] = ys
    preds_by_ratio[ratio] = Dict{String, Vector{Int}}()

    println("  n_total=$(out.stats.n_total) n_low=$(out.stats.n_low) n_high=$(out.stats.n_high) " *
            "high_frac=$(round(100 * out.stats.n_high / out.stats.n_total, digits=1))%")

    for (name, fn) in BACKENDS
        t0 = time()
        pred = fn(smld)
        elapsed_ms = round((time() - t0) * 1000, digits = 1)
        m    = compute_metrics(gt, pred)
        push!(results[name], (ratio = ratio, m...))
        preds_by_ratio[ratio][name] = pred
        @printf("  %-12s acc=%6.2f%% prec=%6.2f%% rec=%6.2f%%   TP=%d FP=%d FN=%d   %.0f ms\n",
                name, 100*m.acc, 100*m.precision, 100*m.recall, m.tp, m.fp, m.fn, elapsed_ms)
    end
end

# ---------------------------------------------------------------------------
# CSV
# ---------------------------------------------------------------------------
csv_path = joinpath(OUT_DIR, "density_ratio_sweep.csv")
open(csv_path, "w") do io
    println(io, "backend,ratio,n_total,accuracy,precision,recall,tp,tn,fp,fn")
    for (name, _) in BACKENDS
        for r in sort(results[name], by = x -> x.ratio)
            n = r.tp + r.tn + r.fp + r.fn
            @printf(io, "%s,%.2f,%d,%.4f,%.4f,%.4f,%d,%d,%d,%d\n",
                    name, r.ratio, n, r.acc, r.precision, r.recall,
                    r.tp, r.tn, r.fp, r.fn)
        end
    end
end
println("[sweep] wrote $csv_path")

# ---------------------------------------------------------------------------
# Aggregate accuracy curve — 4 backends across the 5 ratios.
# ---------------------------------------------------------------------------
fig = Figure(size = (900, 520))
ax  = Axis(fig[1, 1];
           title  = "Density-ratio operational-floor sweep (5×5 μm A431-mimic, fixed patch geometry, kNN k=20)",
           xlabel = "density ratio (high / low)",
           ylabel = "accuracy",
           xticks = (RATIOS, [@sprintf("%.1f×", r) for r in RATIOS]))
hlines!(ax, [0.75]; color = (:gray60, 0.7), linestyle = :dash)
hlines!(ax, [0.85]; color = (:gray40, 0.7), linestyle = :dot)
text!(ax, RATIOS[1] - 0.05, 0.755; text = "75% gate", color = :gray40, fontsize = 11)
text!(ax, RATIOS[1] - 0.05, 0.855; text = "85% gate", color = :gray30, fontsize = 11)
ylims!(ax, 0.4, 1.0)

palette = Dict("mrf_knn"     => :crimson,
               "hdbscan"     => :darkviolet,
               "dbscan"      => :seagreen,
               "voronoi_gmm" => :goldenrod1)

for (name, _) in BACKENDS
    rows  = sort(results[name], by = x -> x.ratio)
    xs_pts = [r.ratio for r in rows]
    accs   = [r.acc   for r in rows]
    lines!(ax, xs_pts, accs;   color = palette[name], linewidth = 2.5, label = name)
    scatter!(ax, xs_pts, accs; color = palette[name], markersize = 12, strokewidth = 0)
end
axislegend(ax; position = :rb, framevisible = false)
save(joinpath(OUT_DIR, "density_ratio_sweep.png"), fig; px_per_unit = 2)
println("[sweep] wrote density_ratio_sweep.png")

# ---------------------------------------------------------------------------
# Categorical-color circle plots — per-emitter, one figure per backend.
# Distinct hues: gray=correct, orange=FP (false high), cyan=FN (missed high).
# ---------------------------------------------------------------------------
const COLOR_OK = (:gray70,    0.45)
const COLOR_FP = (:darkorange, 0.85)
const COLOR_FN = (:darkcyan,   0.85)

function plot_categorical_per_ratio(name)
    fig = Figure(size = (1750, 420))
    Label(fig[0, 1:length(RATIOS)],
          "Per-emitter classification — $name   (gray = correct, orange = FP/false high, cyan = FN/missed high)";
          fontsize = 14)
    for (k, ratio) in enumerate(RATIOS)
        gt   = gt_by_ratio[ratio]
        pred = preds_by_ratio[ratio][name]
        xs   = xs_by_ratio[ratio]
        ys   = ys_by_ratio[ratio]
        m    = compute_metrics(gt, pred)
        ax   = Axis(fig[1, k];
                    title  = @sprintf("%.1f×   acc=%.1f%%   FN=%d  FP=%d",
                                      ratio, 100*m.acc, m.fn, m.fp),
                    xlabel = "x (μm)",
                    ylabel = (k == 1 ? "y (μm)" : ""),
                    aspect = DataAspect())
        xlims!(ax, 0, 5); ylims!(ax, 0, 5)
        colors = [
            gt[i] == pred[i] ? COLOR_OK :
            (gt[i] == 1 && pred[i] == 2 ? COLOR_FP : COLOR_FN)
            for i in eachindex(gt)
        ]
        scatter!(ax, xs, ys; color = colors, markersize = 2.5, strokewidth = 0)
    end
    save(joinpath(OUT_DIR, "density_ratio_sweep_$(name)_panels.png"), fig; px_per_unit = 2)
end

for (name, _) in BACKENDS
    plot_categorical_per_ratio(name)
end

println("\n[sweep] DONE — outputs in $OUT_DIR")
for f in sort(filter(x -> startswith(x, "density_ratio_sweep"), readdir(OUT_DIR)))
    sz = stat(joinpath(OUT_DIR, f)).size
    println("  $f  ($(round(sz/1024, digits=1)) KiB)")
end
