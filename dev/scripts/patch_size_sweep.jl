# Patch-size scaling sweep for the kNN-MRF density-regime clustering pipeline.
# Round 015, Priority 4.
#
# For each patch_scale s ∈ {0.2, 0.33, 0.67, 1.0, 1.33, 2.0} (corresponds to
# mean patch sizes ~{0.3, 0.5, 1.0, 1.5, 2.0, 3.0} μm at the simulator's
# 60% rect / 40% ellipse mix; the default scale=1.0 has rect L ∈ [1,3] μm,
# ellipse semi-major ∈ [0.5,1] μm), regenerate the synthetic A431-mimic at
# that patch scale (same RNG seed so densities/positions are RNG-stable —
# only patch dimensions vary), run kNN-MRF (k=20) and the voronoi-GMM
# baseline, record per-emitter accuracy/precision/recall, and render an
# SMLMRender categorical TP/TN/FP/FN circle plot per scale.
#
# Outputs in dev/scripts/output/:
#   patch_size_sweep.csv                          per (backend, scale) metrics
#   patch_size_sweep.png                          accuracy curve, both backends
#   round_015_smlmrender_categorical.png          6-panel SMLMRender grid (canonical)

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

const OUT_DIR    = joinpath(@__DIR__, "output")
# patch_scale → nominal mean rect length in μm (for labels). simulator's rect L
# range is [scale, 3*scale] so mean ~ 2*scale. We label by 2*scale below.
const PATCH_SCALES = [0.15, 0.25, 0.5, 1.0, 1.5]   # mean rect ~ {0.3, 0.5, 1.0, 2.0, 3.0} μm
const NOMINAL_PATCH_UM = [2 * s for s in PATCH_SCALES]   # ≈ mean rect length

# ---------------------------------------------------------------------------
# Backend predictors (Vector{Int}, 1=low, 2=high) aligned to emitter order.
# ---------------------------------------------------------------------------

function predict_mrf_knn(smld)
    cfg = MRFDensityClusterConfig(n_regimes = 2, density_estimator = :knn, density_k = 20)
    (smld_out, _) = cluster(smld, cfg)
    return smld_out.metadata["mrf_regime_per_emitter"]::Vector{Int}
end

# 2-component GMM EM on log ρ — unsmoothed Voronoi baseline (mirrors Round 013).
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
# Sweep
# ---------------------------------------------------------------------------
const BACKENDS = [
    ("mrf_knn",     predict_mrf_knn),
    ("voronoi_gmm", predict_voronoi_gmm),
]

results = Dict{String, Vector{NamedTuple}}()
for (name, _) in BACKENDS
    results[name] = NamedTuple[]
end

# For each scale, retain SMLD + GT + per-backend prediction so the SMLMRender
# panel grid below can be built without re-running the simulation.
smld_by_scale = Dict{Float64, Any}()
gt_by_scale   = Dict{Float64, Vector{Int}}()
preds_by_scale = Dict{Float64, Dict{String, Vector{Int}}}()
nominal_by_scale = Dict{Float64, Float64}()

for (idx, scale) in enumerate(PATCH_SCALES)
    nominal = NOMINAL_PATCH_UM[idx]
    println("\n[sweep] patch_scale=$scale  (nominal mean rect length ≈ $nominal μm)")
    out = simulate_dataset(patch_scale = scale, verbose = false)
    smld = out.smld; gt = out.ground_truth
    nominal_by_scale[scale] = nominal
    smld_by_scale[scale]    = smld
    gt_by_scale[scale]      = gt
    preds_by_scale[scale]   = Dict{String, Vector{Int}}()

    println("  n_total=$(out.stats.n_total) n_low=$(out.stats.n_low) n_high=$(out.stats.n_high) " *
            "high_frac=$(round(100 * out.stats.n_high / out.stats.n_total, digits=1))%  " *
            "n_patches=$(out.stats.n_patches)")

    for (name, fn) in BACKENDS
        t0 = time()
        pred = fn(smld)
        elapsed_ms = round((time() - t0) * 1000, digits = 1)
        m = compute_metrics(gt, pred)
        push!(results[name], (scale = scale, nominal_um = nominal, m...))
        preds_by_scale[scale][name] = pred
        @printf("  %-12s acc=%6.2f%% prec=%6.2f%% rec=%6.2f%%   TP=%d FP=%d FN=%d   %.0f ms\n",
                name, 100*m.acc, 100*m.precision, 100*m.recall, m.tp, m.fp, m.fn, elapsed_ms)
    end
end

# ---------------------------------------------------------------------------
# CSV
# ---------------------------------------------------------------------------
csv_path = joinpath(OUT_DIR, "patch_size_sweep.csv")
open(csv_path, "w") do io
    println(io, "backend,patch_scale,nominal_um,n_total,accuracy,precision,recall,tp,tn,fp,fn")
    for (name, _) in BACKENDS
        for r in sort(results[name], by = x -> x.scale)
            n = r.tp + r.tn + r.fp + r.fn
            @printf(io, "%s,%.3f,%.3f,%d,%.4f,%.4f,%.4f,%d,%d,%d,%d\n",
                    name, r.scale, r.nominal_um, n, r.acc, r.precision, r.recall,
                    r.tp, r.tn, r.fp, r.fn)
        end
    end
end
println("\n[sweep] wrote $csv_path")

# ---------------------------------------------------------------------------
# Aggregate accuracy curve — both backends across the patch-size sweep.
# ---------------------------------------------------------------------------
fig = Figure(size = (900, 520))
ax  = Axis(fig[1, 1];
           title  = "Patch-size scaling sweep (5×5 μm A431-mimic, 2× density ratio, kNN k=20)",
           xlabel = "nominal patch size (μm — mean rect length, log scale)",
           ylabel = "accuracy",
           xscale = log10,
           xticks = (NOMINAL_PATCH_UM, [@sprintf("%.1f", n) for n in NOMINAL_PATCH_UM]))
hlines!(ax, [0.75]; color = (:gray60, 0.7), linestyle = :dash)
hlines!(ax, [0.85]; color = (:gray40, 0.7), linestyle = :dot)
text!(ax, NOMINAL_PATCH_UM[1] * 0.85, 0.755; text = "75% gate", color = :gray40, fontsize = 11)
text!(ax, NOMINAL_PATCH_UM[1] * 0.85, 0.855; text = "85% gate", color = :gray30, fontsize = 11)
ylims!(ax, 0.4, 1.0)

palette = Dict("mrf_knn"     => :crimson,
               "voronoi_gmm" => :goldenrod1)

for (name, _) in BACKENDS
    rows = sort(results[name], by = x -> x.scale)
    xs_pts = [r.nominal_um for r in rows]
    accs   = [r.acc for r in rows]
    lines!(ax, xs_pts, accs;   color = palette[name], linewidth = 2.5, label = name)
    scatter!(ax, xs_pts, accs; color = palette[name], markersize = 12, strokewidth = 0)
end
axislegend(ax; position = :rb, framevisible = false)
save(joinpath(OUT_DIR, "patch_size_sweep.png"), fig; px_per_unit = 2)
println("[sweep] wrote patch_size_sweep.png")

# ---------------------------------------------------------------------------
# SMLMRender categorical-color circle plot — canonical Round 015 deliverable.
# 5 panels, one per patch-size scale, showing kNN-MRF TP/TN/FP/FN per emitter.
# Uses SMLMRender.render(..., strategy=CircleRender(), color_by=:id, categorical=true).
# ---------------------------------------------------------------------------

# Encode TP/TN/FP/FN as integer ids 1..4 by stamping a deepcopy SMLD's emitter.id:
#   1 = TN (correct low)        gray
#   2 = TP (correct high)       green
#   3 = FP (predicted high but actually low)  orange
#   4 = FN (predicted low but actually high)  cyan
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

# tab10 first four hues (categorical, qualitative) — but we want
# deterministic semantic colors, so use ManualColorMapping. SMLMRender's
# CategoricalColorMapping with :tab10 cycles ids 1..N to palette[mod1(id, N)].
# For 4 distinct semantic categories, :tab10 gives blue/orange/green/red —
# good enough; we relabel ids so the semantic-to-color mapping is clear.
# id_to_tab10:  1 (TN) → tab10[1] blue  → relabel TN ids to a less salient slot
# Simpler: use :Set1_9 palette — first 4 hues are red, blue, green, purple.
# Even simpler: re-map our ids so the qualitative palette gives the right vibe.
# Convention: stamp ids 1=TN(faint), 2=TP, 3=FP, 4=FN. CategoricalColorMapping
# uses :tab10 default — distinct hues per id is the requirement.

function render_per_scale_panel(smld, gt, pred; zoom::Int = 8)
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
panel_imgs = Array{Any}(undef, length(PATCH_SCALES))
for (idx, scale) in enumerate(PATCH_SCALES)
    smld = smld_by_scale[scale]
    gt   = gt_by_scale[scale]
    pred = preds_by_scale[scale]["mrf_knn"]
    panel_imgs[idx] = render_per_scale_panel(smld, gt, pred; zoom = 8)
    println("  scale=$scale  rendered  size=$(size(panel_imgs[idx]))")
end

# Compose into a single 1xN multipanel figure via CairoMakie, image! per axis.
fig2 = Figure(size = (1850, 460))
Label(fig2[0, 1:length(PATCH_SCALES)],
      "kNN-MRF (k=20) per-emitter classification — patch-size sweep   "
      * "(SMLMRender CircleRender, categorical: id 1=TN, 2=TP, 3=FP, 4=FN)";
      fontsize = 14)
for (k, scale) in enumerate(PATCH_SCALES)
    nominal = nominal_by_scale[scale]
    gt   = gt_by_scale[scale]
    pred = preds_by_scale[scale]["mrf_knn"]
    m = compute_metrics(gt, pred)
    img = panel_imgs[k]
    ax = Axis(fig2[1, k];
              title = @sprintf("~%.1f μm   acc=%.1f%%   FN=%d  FP=%d",
                               nominal, 100*m.acc, m.fn, m.fp),
              aspect = DataAspect())
    hidedecorations!(ax)
    image!(ax, rotr90(img))
end
save(joinpath(OUT_DIR, "round_015_smlmrender_categorical.png"), fig2; px_per_unit = 2)
println("[smlmrender] wrote round_015_smlmrender_categorical.png")

println("\n[sweep] DONE — outputs in $OUT_DIR")
for f in sort(filter(x -> startswith(x, "patch_size_sweep") || startswith(x, "round_015"),
                     readdir(OUT_DIR)))
    sz = stat(joinpath(OUT_DIR, f)).size
    println("  $f  ($(round(sz/1024, digits=1)) KiB)")
end
