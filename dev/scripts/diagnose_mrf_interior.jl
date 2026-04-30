# Diagnose the MRF interior false-negative band on the synthetic A431-mimic.
#
# Round 012 / Priority 1.
#
# Hypothesis: MRF v1 (Voronoi density) produces a concentrated cyan FN band
# at high-density patch interiors because the Voronoi cell linear scale at
# 1000/μm² (~30-40 nm) is comparable to thin-fiber widths (50-600 nm for AR
# 5-20 rectangles). Cells that span the boundary inflate, dragging interior
# log-densities toward the low regime.
#
# Mitigation under test: kNN density estimator (`density_estimator=:knn` with
# `density_k=20`), which integrates over the k nearest neighbors and reduces
# noise from σ_log≈1.0 (single-cell Voronoi) to σ_log≈1/√k≈0.22 (k=20).
#
# Outputs:
#   dev/scripts/output/mrf_interior_diagnosis.png      # 6-panel figure
#   dev/scripts/output/mrf_interior_diagnosis.csv      # per-patch metrics

using Pkg
Pkg.activate(@__DIR__)

using SMLMClustering
using SMLMData
using JLD2
using CairoMakie
using NearestNeighbors
using Statistics

const OUT_DIR = joinpath(@__DIR__, "output")
const ROI_X_HI = 4.375
const ROI_Y_HI = 5.0
const KNN_K = 20

# ---------------------------------------------------------------------------
# Load
# ---------------------------------------------------------------------------
println("[diag] loading synthetic SMLD…")
data = load(joinpath(OUT_DIR, "synthetic_smld.jld2"))
smld     = data["smld"]
gt       = data["ground_truth"]::Vector{Int}
patches  = data["patches"]
n        = length(smld.emitters)
xs       = [e.x for e in smld.emitters]
ys       = [e.y for e in smld.emitters]
println("  n=$n | n_low=$(count(==(1), gt)) | n_high=$(count(==(2), gt)) | n_patches=$(length(patches))")

# ---------------------------------------------------------------------------
# Per-emitter density estimates.
# ---------------------------------------------------------------------------
println("[diag] computing per-emitter densities (Voronoi via cluster_statistics + kNN k=$KNN_K)…")
(_, vinfo) = cluster_statistics(smld, VoronoiDensityConfig())
ρ_voronoi   = vinfo.extras[:density_per_emitter]::Vector{Float64}

# kNN density: ρ_k = k / (π · r_k²); replicates SMLMClustering._knn_density.
function knn_density(xs::Vector{Float64}, ys::Vector{Float64}, k::Int)
    n = length(xs)
    X = Matrix{Float64}(undef, 2, n)
    @inbounds for i in 1:n; X[1,i]=xs[i]; X[2,i]=ys[i]; end
    tree = KDTree(X)
    out = Vector{Float64}(undef, n)
    @inbounds for i in 1:n
        idxs, dists = knn(tree, view(X, :, i), k+1, true)
        rk = dists[end]
        out[i] = rk > 0 ? k / (π * rk^2) : NaN
    end
    return out
end
ρ_knn = knn_density(xs, ys, KNN_K)

# ---------------------------------------------------------------------------
# Run MRF: :voronoi (current default) and :knn (mitigation under test).
# ---------------------------------------------------------------------------
println("[diag] running MRF :voronoi…")
cfg_v = MRFDensityClusterConfig(n_regimes=2, density_estimator=:voronoi)
(out_v, _) = cluster(smld, cfg_v)
regime_v = out_v.metadata["mrf_regime_per_emitter"]::Vector{Int}

println("[diag] running MRF :knn (k=$KNN_K)…")
cfg_k = MRFDensityClusterConfig(n_regimes=2, density_estimator=:knn, density_k=KNN_K)
(out_k, _) = cluster(smld, cfg_k)
regime_k = out_k.metadata["mrf_regime_per_emitter"]::Vector{Int}

# Convert regime → 2-class label (1=low, 2=high). Regime 0 (ungroupable) → 1.
predict_label(r) = r == 2 ? 2 : 1
pred_v = predict_label.(regime_v)
pred_k = predict_label.(regime_k)

# ---------------------------------------------------------------------------
# Confusion + accuracy
# ---------------------------------------------------------------------------
function confusion(gt, pred)
    a=b=c=d=0
    for i in eachindex(gt)
        if gt[i]==1 && pred[i]==1; a+=1
        elseif gt[i]==1 && pred[i]==2; b+=1
        elseif gt[i]==2 && pred[i]==1; c+=1
        else; d+=1; end
    end
    (low_low=a, low_high=b, high_low=c, high_high=d, acc=(a+d)/length(gt))
end
cv = confusion(gt, pred_v)
ck = confusion(gt, pred_k)

println("\n[diag] HEADLINE ACCURACY")
println("  voronoi : ", round(100*cv.acc, digits=2), "%   FP=", cv.low_high, "  FN=", cv.high_low)
println("  knn(k=$KNN_K): ", round(100*ck.acc, digits=2), "%   FP=", ck.low_high, "  FN=", ck.high_low)

# ---------------------------------------------------------------------------
# Distance-to-patch-boundary (signed: negative inside patches, positive outside).
# Approximate using the patch geometry primitives. For a point at (x,y) and a
# patch p, the in_patch test gives sign; magnitude = closest perpendicular
# distance in the patch-local frame. Diagnostic-grade approximation; for
# rectangles uses dist-to-rect; for ellipses uses normalized-radius approx.
# ---------------------------------------------------------------------------
function patch_local(p, x, y)
    dx = x - p.cx; dy = y - p.cy
    c = cos(-p.theta); s = sin(-p.theta)
    return (c*dx - s*dy, s*dx + c*dy)
end

function signed_dist_to_patch(p, x, y)
    xl, yl = patch_local(p, x, y)
    if p.kind === :rect
        # distance to rectangle [-a,a]×[-b,b]: positive outside, negative inside.
        ax = abs(xl) - p.a; ay = abs(yl) - p.b
        if ax <= 0 && ay <= 0
            # interior — distance to nearest side, sign negative
            return -min(-ax, -ay)
        else
            ax_p = max(ax, 0.0); ay_p = max(ay, 0.0)
            return sqrt(ax_p^2 + ay_p^2)
        end
    else  # :ellipse
        # Use normalized radius to label inside/outside; magnitude is approximate.
        nr = sqrt((xl/p.a)^2 + (yl/p.b)^2)
        # Linear approximation: dist ≈ (nr - 1) * min(p.a, p.b)
        return (nr - 1) * min(p.a, p.b)
    end
end

# Per-emitter signed distance to NEAREST patch (negative = interior).
println("[diag] computing per-emitter distance to nearest patch boundary…")
sd = Vector{Float64}(undef, n)
nearest_patch = zeros(Int, n)
@inbounds for i in 1:n
    best_d = Inf
    best_p = 0
    best_signed = Inf
    inside_any = false
    for (k, p) in enumerate(patches)
        d = signed_dist_to_patch(p, xs[i], ys[i])
        # If inside this patch, take the "most-interior" (most-negative) one.
        if d < 0
            inside_any = true
            if d < best_signed
                best_signed = d
                best_p = k
            end
        elseif !inside_any && d < best_d
            best_d = d
            best_p = k
            best_signed = d
        end
    end
    sd[i] = best_signed
    nearest_patch[i] = best_p
end

# Interior (sd < -100 nm = -0.1 μm) is the diagnosis target.
interior_mask = sd .< -0.1
println("  interior emitters (within -100 nm of nearest patch boundary): ",
        count(interior_mask), " of ", n, " (",
        round(100*count(interior_mask)/n, digits=1), "%)")

# Per-region FN-rate
function fn_rate(gt, pred, mask)
    g_high_in_mask = (gt .== 2) .& mask
    fn_in_mask = (gt .== 2) .& (pred .== 1) .& mask
    n_high = count(g_high_in_mask)
    n_high == 0 ? NaN : count(fn_in_mask) / n_high
end

fn_int_v = fn_rate(gt, pred_v, interior_mask)
fn_int_k = fn_rate(gt, pred_k, interior_mask)
fn_all_v = fn_rate(gt, pred_v, trues(n))
fn_all_k = fn_rate(gt, pred_k, trues(n))

println("\n[diag] PATCH-INTERIOR FN-RATE  (high-density emitters >100 nm inside a patch)")
println("  voronoi  interior FN-rate: ", round(100*fn_int_v, digits=2), "%   (overall ",
        round(100*fn_all_v, digits=2), "%)")
println("  knn      interior FN-rate: ", round(100*fn_int_k, digits=2), "%   (overall ",
        round(100*fn_all_k, digits=2), "%)")

# ---------------------------------------------------------------------------
# Per-patch metrics → CSV
# ---------------------------------------------------------------------------
println("[diag] writing per-patch metrics CSV…")
csv_path = joinpath(OUT_DIR, "mrf_interior_diagnosis.csv")
open(csv_path, "w") do io
    println(io, "patch,kind,a_um,b_um,aspect_ratio,area_um2,n_high_interior,",
            "fn_int_voronoi,fn_int_knn")
    for (k, p) in enumerate(patches)
        pmask = (nearest_patch .== k) .& interior_mask
        n_hi = count((gt .== 2) .& pmask)
        if n_hi == 0
            v_rate = NaN; k_rate = NaN
        else
            v_rate = count((gt .== 2) .& (pred_v .== 1) .& pmask) / n_hi
            k_rate = count((gt .== 2) .& (pred_k .== 1) .& pmask) / n_hi
        end
        ar = p.a / max(p.b, 1e-12)
        println(io, k, ",", p.kind, ",", round(p.a; digits=4), ",",
                round(p.b; digits=4), ",", round(ar; digits=2), ",",
                round(p.area; digits=4), ",", n_hi, ",",
                round(v_rate; digits=4), ",", round(k_rate; digits=4))
    end
end
println("  saved $csv_path")

# ---------------------------------------------------------------------------
# Plotting helpers
# ---------------------------------------------------------------------------
const COLOR_LOW   = (:steelblue, 0.6)
const COLOR_HIGH  = (:crimson,   0.6)
const COLOR_OK    = (:gray70,    0.4)
const COLOR_FP    = (:darkorange, 0.8)
const COLOR_FN    = (:darkcyan,   0.8)

function patch_outline_pts(p)
    if p.kind === :rect
        c=cos(p.theta); s=sin(p.theta)
        loc = [(p.a, p.b), (-p.a, p.b), (-p.a, -p.b), (p.a, -p.b), (p.a, p.b)]
        return [(p.cx + c*x - s*y, p.cy + s*x + c*y) for (x,y) in loc]
    else
        c=cos(p.theta); s=sin(p.theta)
        return [
            begin
                t = 2π * (kk - 1) / 64
                xl = p.a*cos(t); yl = p.b*sin(t)
                (p.cx + c*xl - s*yl, p.cy + s*xl + c*yl)
            end
            for kk in 1:65
        ]
    end
end

function draw_patches!(ax, patches)
    for p in patches
        pts = patch_outline_pts(p)
        lines!(ax, [t[1] for t in pts], [t[2] for t in pts];
               color=:black, linewidth=0.5)
    end
end

function shade_edge!(ax)
    poly!(ax, [(ROI_X_HI, 0.0), (ROI_Y_HI, 0.0),
               (ROI_Y_HI, ROI_Y_HI), (ROI_X_HI, ROI_Y_HI)];
          color=(:gray80, 0.3), strokewidth=0)
end

# ---------------------------------------------------------------------------
# Main 6-panel figure
#   row 1: log-density distributions (Voronoi vs kNN), colored by GT
#   row 2: error maps (Voronoi vs kNN); cyan=FN, orange=FP, gray=correct
#   row 3: signed-distance histogram of FN events (Voronoi vs kNN)
# ---------------------------------------------------------------------------
println("[diag] rendering 6-panel diagnostic figure…")
fig = Figure(size=(1500, 1500))

# --- Row 1: density distributions ---
log_ρ_v = log.(filter(isfinite, ρ_voronoi))
log_ρ_k = log.(filter(isfinite, ρ_knn))
gt_v = gt[isfinite.(ρ_voronoi)]
gt_k = gt[isfinite.(ρ_knn)]

ax1a = Axis(fig[1,1]; title="Voronoi log-density (Voronoi 1/A)",
            xlabel="log ρ (μm⁻²)", ylabel="count")
hist!(ax1a, log_ρ_v[gt_v.==1]; bins=80, color=(:steelblue, 0.5), label="GT low")
hist!(ax1a, log_ρ_v[gt_v.==2]; bins=80, color=(:crimson, 0.5),  label="GT high")
axislegend(ax1a; position=:lt)

ax1b = Axis(fig[1,2]; title="kNN log-density (k=$KNN_K)",
            xlabel="log ρ (μm⁻²)", ylabel="count")
hist!(ax1b, log_ρ_k[gt_k.==1]; bins=80, color=(:steelblue, 0.5), label="GT low")
hist!(ax1b, log_ρ_k[gt_k.==2]; bins=80, color=(:crimson, 0.5),  label="GT high")
axislegend(ax1b; position=:lt)

# --- Row 2: error maps ---
function err_color(gti, predi)
    gti == predi && return COLOR_OK
    (gti == 1 && predi == 2) ? COLOR_FP : COLOR_FN
end

ax2a = Axis(fig[2,1];
            title="MRF :voronoi errors (acc=$(round(100*cv.acc, digits=1))%, FN=$(cv.high_low))",
            xlabel="x (μm)", ylabel="y (μm)", aspect=DataAspect())
xlims!(ax2a, 0, 5); ylims!(ax2a, 0, 5)
scatter!(ax2a, xs, ys; color=[err_color(gt[i], pred_v[i]) for i in 1:n],
         markersize=2.5, strokewidth=0)
draw_patches!(ax2a, patches); shade_edge!(ax2a)

ax2b = Axis(fig[2,2];
            title="MRF :knn errors (acc=$(round(100*ck.acc, digits=1))%, FN=$(ck.high_low))",
            xlabel="x (μm)", ylabel="y (μm)", aspect=DataAspect())
xlims!(ax2b, 0, 5); ylims!(ax2b, 0, 5)
scatter!(ax2b, xs, ys; color=[err_color(gt[i], pred_k[i]) for i in 1:n],
         markersize=2.5, strokewidth=0)
draw_patches!(ax2b, patches); shade_edge!(ax2b)

# --- Row 3: signed-distance histograms of FN events ---
fn_mask_v = (gt .== 2) .& (pred_v .== 1)
fn_mask_k = (gt .== 2) .& (pred_k .== 1)
ax3a = Axis(fig[3,1];
            title="Voronoi: FN locations vs distance from patch boundary",
            xlabel="signed distance to nearest patch boundary (μm)\n(negative = inside patch)",
            ylabel="count of FN events")
hist!(ax3a, sd[fn_mask_v]; bins=40, color=(:darkcyan, 0.7))
vlines!(ax3a, 0; color=:black, linestyle=:dash)
vlines!(ax3a, -0.1; color=:gray50, linestyle=:dot)
text!(ax3a, -0.1, 1; text=" -100 nm threshold", fontsize=10, color=:gray50)

ax3b = Axis(fig[3,2];
            title="kNN: FN locations vs distance from patch boundary",
            xlabel="signed distance to nearest patch boundary (μm)\n(negative = inside patch)",
            ylabel="count of FN events")
hist!(ax3b, sd[fn_mask_k]; bins=40, color=(:darkcyan, 0.7))
vlines!(ax3b, 0; color=:black, linestyle=:dash)
vlines!(ax3b, -0.1; color=:gray50, linestyle=:dot)
text!(ax3b, -0.1, 1; text=" -100 nm threshold", fontsize=10, color=:gray50)

Label(fig[0, :], "MRF interior FN diagnosis — synthetic A431-mimic | Round 012 P1";
      fontsize=18, font=:bold)

png_path = joinpath(OUT_DIR, "mrf_interior_diagnosis.png")
save(png_path, fig; px_per_unit=2)
println("[diag] saved $png_path")

# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------
println("\n=== SUMMARY ===")
println("Voronoi:  acc=", round(100*cv.acc, digits=2), "%   interior FN-rate=",
        round(100*fn_int_v, digits=2), "%")
println("kNN(k=$KNN_K): acc=", round(100*ck.acc, digits=2), "%   interior FN-rate=",
        round(100*fn_int_k, digits=2), "%")
println("Δ headline accuracy: +", round(100*(ck.acc-cv.acc), digits=2), " pp")
println("Δ interior FN-rate:  ", round(100*(fn_int_k-fn_int_v), digits=2), " pp")

ck.acc > 0.85 || @warn "kNN headline accuracy below 85% — round may close as `partial`"
fn_int_k < fn_int_v ||
    @warn "kNN did NOT reduce interior FN-rate — hypothesis inverted; needs second look"
