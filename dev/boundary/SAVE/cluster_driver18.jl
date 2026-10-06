using SMLMAnalysis        # load_smld
using SMLMData            # filter_roi
using SMLMClustering      # cluster, cluster_statistics, boundary_clusters, …
using Statistics          # mean, std

include("plot_clusters.jl")   # ClusterPlotConfig, plot_clusters

## ─────────────────────────────────────────────────────────────────────────────
## Main driver
## ─────────────────────────────────────────────────────────────────────────────

# Path to input h5 smld dataset
#smld_path = "/mnt/nas/cellpath/papers-geometry/julia-analysis-results/dstorm/20231218_IgE-DF3/DF310nM5min/Cell_06/14_render/gaussianrender_tab10_id_50x.png"
smld_path = "/mnt/nas/lidkelab/Personal Folders/MJW/Julia/publish/resultsH5/" *
            "Data_2025-11-20-10-52-19_result_smld.h5"
smld = load_smld(smld_path)
smld = filter_roi(smld, (12.5, 17.5), (12.5, 17.5))   # 5 μm × 5 μm centered at (15, 15) μm
#smld = filter_roi(smld, (12.5, 13.5), (12.5, 13.5))

# Path to directory to place results.
results_dir = "/mnt/nas/lidkelab/Personal Folders/MJW/Julia/"

# Cluster within each dataset independently; must match in the configuration
# parameters for cluster and boundary_clusters,
per_dataset = false

# rendered pixel size (nm); smaller = finer detail
pixel_size_nm = 1.0

# ── Clustering (DBSCAN) ───────────────────────────────────────────────────────
cfg_dbscan = DBSCANConfig(
    eps_nm             = 25.0,  # neighborhood radius in nm (required)
    min_points         = 3,     # core-point threshold / min cluster size
    use_3d             = false, # include z-coordinate
    per_dataset        = per_dataset, # cluster within each dataset independently
    remove_unclustered = false,
)
(smld_out, info_dbscan) = cluster(smld, cfg_dbscan)
println(smld_out)
println(info_dbscan)

# ── Cluster boundaries (alpha-shape) ─────────────────────────────────────────
cfg_bnd = BoundaryClustersConfig(
    shrink      = 0.05,  # 0 = convex hull, 1 = tightest alpha-shape
    per_dataset = per_dataset, # must match the per_dataset value used in cluster()
)
(_, binfo) = boundary_clusters(smld_out, cfg_bnd)
println(binfo)
# binfo.boundaries_x[j] / binfo.boundaries_y[j] trace the boundary polygon of
# cluster binfo.cluster_keys[j] in μm (closed: first point == last point).

# ── Cluster boundary statistics ───────────────────────────────────────────────
(_, sinfo) = boundary_statistics(smld_out, binfo)
println(sinfo)

valid_areas  = filter(!isnan, sinfo.areas)
valid_perims = filter(!isnan, sinfo.perimeters)

if !isempty(valid_areas)
    n_a = length(valid_areas)
    println("Area:      = ", round(mean(valid_areas), digits = 4),
            "  +/- ", round(std(valid_areas) / sqrt(n_a), digits = 4), " μm²",
            "  (n = ", n_a, ")")
end
if !isempty(valid_perims)
    n_p = length(valid_perims)
    println("Perimeter: = ", round(mean(valid_perims), digits = 4),
            "  +/- ", round(std(valid_perims) / sqrt(n_p), digits = 4), " μm ",
            "  (n = ", n_p, ")")
end

# ── Cluster visualization ─────────────────────────────────────────────────────
# Derive a PNG filename from the data file (placed next to this script).
png_path = results_dir * splitext(basename(smld_path))[1] * "_Gaussian_clusters.png"
plot_cfg = ClusterPlotConfig(
    pixel_size_nm = pixel_size_nm,
    colormap      = :inferno,      # good contrast for SMLM density images
    color_by      = :id,           # color each localization by its cluster id
    categorical   = true,
    boundaries    = binfo,         # draw alpha-shape outlines for each cluster
    show_centroids = false,
    filename      = png_path,
)
t_plot = @elapsed fig = plot_clusters(smld, smld_out, plot_cfg)
println("plot_clusters: ", round(t_plot, digits = 3), " s")

# ── Cluster visualization (fixed-sigma "dot" render) ──────────────────────────
png_path2 = results_dir * splitext(basename(smld_path))[1] * "_dot_clusters.png"
plot_cfg2 = ClusterPlotConfig(
    pixel_size_nm   = pixel_size_nm,
    colormap        = :inferno,
    render_strategy = GaussianRender(use_localization_precision = false,
                                      fixed_sigma = 1.0),
    boundaries      = binfo,
    boundary_color  = :white,
    show_centroids  = false,
    filename        = png_path2,
)
t_plot2 = @elapsed fig2 = plot_clusters(smld, smld_out, plot_cfg2)
println("plot_clusters (dot): ", round(t_plot2, digits = 3), " s")
