# Voronoi-tessellation backend (SR-Tesseler-style density clustering).
#
# Localizations are clustered by local Voronoi-cell density: each point's cell
# area is compared against the group's mean area; points in cells smaller than
# `mean_area / density_factor` are "dense" and get agglomerated into clusters
# via connected components over the Delaunay adjacency graph. Small components
# are relabeled noise via `min_points`.
#
# Reference: Levet et al., "SR-Tesseler: a method to segment and quantify
# localization-based super-resolution microscopy data," Nat. Methods (2015).
#
# Uses DelaunayTriangulation.jl (pure Julia, 2D). 3D clustering is not
# supported — `use_3d=true` raises an error.

using DelaunayTriangulation

"""
    VoronoiConfig(; density_factor=2.0, min_points=5, use_3d=false,
                    per_dataset=true, remove_unclustered=false)

Configuration for Voronoi-tessellation-based (SR-Tesseler) clustering of SMLM
localizations.

# Algorithm
1. Build the Voronoi tessellation of the localization coordinates (per
   dataset when `per_dataset=true`), clipped to the convex hull so every
   generator has a finite cell.
2. A localization is **dense** when its cell area is smaller than
   `mean_cell_area / density_factor` (equivalently, its local density exceeds
   `density_factor × mean_density`).
3. Dense points that are Delaunay-adjacent are merged into clusters via
   connected components.
4. Clusters with fewer than `min_points` members are relabeled noise
   (`id = 0`); the rest are renumbered compactly `1..K` within the group.

# Fields
- `density_factor::Float64 = 2.0`: density threshold multiplier. Higher values
  require stronger local density to qualify as a cluster member (smaller area
  threshold = fewer dense points).
- `min_points::Int = 5`: minimum cluster size; smaller connected components
  become noise.
- `use_3d::Bool = false`: must be `false`. 3D Voronoi clustering is not
  supported — DelaunayTriangulation.jl is 2D only.
- `per_dataset::Bool = true`: cluster within each dataset independently so
  that `(dataset, id)` uniquely identifies a cluster.
- `remove_unclustered::Bool = false`: drop noise emitters (`id == 0`) from the
  returned SMLD.

!!! note "Boundary handling"
    Cells are clipped to the convex hull of the generator set. Generators on
    the hull get cells smaller than their true infinite-plane area, which can
    bias mean-area estimates on very small groups. In practice the effect is
    second-order for SMLM datasets with thousands of localizations.

!!! note "Degenerate input"
    Groups with fewer than 3 points are tagged all-noise (a tessellation
    requires at least 3 non-collinear points). Groups containing exact-duplicate
    (x,y) coordinate pairs raise `ArgumentError`; deduplicate input
    localizations before calling `cluster`.

# Example
```julia
cfg = VoronoiConfig(density_factor=2.0, min_points=5)
(smld_out, info) = cluster(smld, cfg)
```

See also: [`AbstractClusterConfig`](@ref), [`ClusterInfo`](@ref), [`cluster`](@ref).
"""
Base.@kwdef struct VoronoiConfig <: AbstractClusterConfig
    density_factor::Float64 = 2.0
    min_points::Int = 5
    use_3d::Bool = false
    per_dataset::Bool = true
    remove_unclustered::Bool = false
end

function cluster(smld::SMLMData.BasicSMLD, cfg::VoronoiConfig)
    t0 = time_ns()
    # Non-mutating semantics: deep-copy emitters so cluster labels go to a
    # fresh SMLD, not back onto the caller's input. See KB V9.
    smld = SMLMData.BasicSMLD(deepcopy(smld.emitters), smld.camera,
                              smld.n_frames, smld.n_datasets, smld.metadata)
    n_in = length(smld.emitters)
    cfg.density_factor > 0 ||
        throw(ArgumentError("VoronoiConfig.density_factor must be > 0 (got $(cfg.density_factor))"))
    cfg.min_points >= 1 ||
        throw(ArgumentError("VoronoiConfig.min_points must be ≥ 1 (got $(cfg.min_points))"))
    cfg.use_3d &&
        throw(ArgumentError(
            "VoronoiConfig does not support use_3d=true. " *
            "3D Voronoi tessellation is not implemented by DelaunayTriangulation.jl; " *
            "use DBSCANConfig or HierarchicalConfig with use_3d=true for 3D data."))

    groups = _group_by_dataset(smld, cfg.per_dataset)

    cluster_sizes = Int[]
    n_clustered = 0

    for idxs in groups
        n = length(idxs)
        if n < 3
            # Tessellation requires ≥3 non-collinear points; tag all as noise.
            @inbounds for i in idxs
                smld.emitters[i].id = 0
            end
            continue
        end

        sub = view(smld.emitters, idxs)
        pts = [(sub[j].x, sub[j].y) for j in 1:n]  # μm

        # Exact-coincident generators cause get_area to raise KeyError.
        length(unique(pts)) == n ||
            throw(ArgumentError(
                "VoronoiConfig: group of $n points contains duplicate (x,y) " *
                "coordinates; deduplicate input localizations before calling cluster()."))

        tri = DelaunayTriangulation.triangulate(pts)
        vor = DelaunayTriangulation.voronoi(tri; clip = true)

        areas = Vector{Float64}(undef, n)
        @inbounds for j in 1:n
            areas[j] = DelaunayTriangulation.get_area(vor, j)
        end

        # Dense ⇔ cell area < mean_area / density_factor.
        mean_area = sum(areas) / n
        area_thresh = mean_area / cfg.density_factor
        dense = Vector{Bool}(undef, n)
        @inbounds for j in 1:n
            dense[j] = areas[j] < area_thresh
        end

        # Connected components on dense points via Delaunay adjacency
        # (ghost neighbour -1 is filtered).
        raw_labels = zeros(Int, n)
        k_raw = 0
        stack = Int[]
        @inbounds for seed in 1:n
            (dense[seed] && raw_labels[seed] == 0) || continue
            k_raw += 1
            raw_labels[seed] = k_raw
            push!(stack, seed)
            while !isempty(stack)
                v = pop!(stack)
                for w in DelaunayTriangulation.get_neighbours(tri, v)
                    if w > 0 && dense[w] && raw_labels[w] == 0
                        raw_labels[w] = k_raw
                        push!(stack, w)
                    end
                end
            end
        end

        # Tally raw component sizes, drop those below min_points, compact relabel.
        raw_counts = zeros(Int, k_raw)
        @inbounds for l in raw_labels
            l > 0 && (raw_counts[l] += 1)
        end
        label_map, added = _compact_relabel!(cluster_sizes, raw_counts, cfg.min_points)
        n_clustered += added

        # Write final labels (local to group per V3).
        @inbounds for (j, i) in pairs(idxs)
            rl = raw_labels[j]
            smld.emitters[i].id = rl > 0 ? label_map[rl] : 0
        end
    end

    n_clusters = length(cluster_sizes)
    n_noise = n_in - n_clustered
    smld_out = _build_output(smld, cfg.remove_unclustered)

    info = ClusterInfo(
        n_in,
        n_clustered,
        n_noise,
        n_clusters,
        cluster_sizes,
        :voronoi,
        (time_ns() - t0) / 1e9,
    )
    return smld_out, info
end
