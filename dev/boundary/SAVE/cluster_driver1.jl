using SMLMAnalysis
using SMLMClustering

smld_path = "/mnt/nas/lidkelab/Personal Folders/MJW/Julia/publish/resultsH5/Data_2025-11-20-10-52-19_result_smld.h5"
smld = load_smld(smld_path)

cfg = DBSCANConfig(
    eps_nm          = 50.0,   # neighborhood radius in nm (required)
    min_points      = 5,      # core-point threshold / min cluster size
    use_3d          = false,  # include z-coordinate
    per_dataset     = true,   # cluster within each dataset independently
    remove_unclustered = false,
)
(smld_out, info) = cluster(smld, cfg)
println(smld_out)
println(info)

cfg_bnd = BoundaryClustersConfig(
    shrink      = 0.5,   # 0 = convex hull, 1 = tightest alpha-shape
    per_dataset = true,  # match the per_dataset value used in cluster()
)
(_, binfo) = boundary_clusters(smld_out, cfg_bnd)
println(binfo)
# binfo.boundaries_x[j], binfo.boundaries_y[j] trace the boundary polygon of
# cluster binfo.cluster_keys[j] in μm (closed: first point == last point).

## Save a PNG next to this script; derive the name from the data filename
#png_path = splitext(basename(smld_path))[1] * "_dbscan.png"
#save(png_path, fig)
#println("Saved → $png_path")

cfg = HopkinsConfig(
    n_samples       = 50,        # reference / sampled point count per repeat
    random_repeats  = 10,        # average over independent repeats
    seed            = 1,         # RNG seed for reproducibility
    use_3d          = false,
    per_dataset     = true,      # report per-dataset H in extras + mean as `statistic`
)
(_, info) = cluster_statistics(smld, cfg)
println("Hopkins H = ", round(info.statistic, digits = 3))
println("per-dataset H = ", info.extras[:hopkins_per_dataset])
