# HDBSCAN backend.
#
# Hierarchical Density-Based Spatial Clustering of Applications with Noise
# (Campello, Moulavi, Sander 2013). Pure-Julia implementation built on top
# of `NearestNeighbors.KDTree`. Configuration subtypes `AbstractClusterConfig`;
# `cluster(smld, ::HDBSCANConfig)` writes per-emitter cluster labels into
# `emitter.id` (0 = noise, 1..K = cluster) and returns `(smld_out, ClusterInfo)`.
#
# HDBSCAN-specific outputs (per-cluster persistence/stability and birth lambda)
# live in `smld_out.metadata`:
#   - `metadata["hdbscan_cluster_persistence"]::Vector{Float64}` (length K)
#   - `metadata["hdbscan_cluster_lambda_birth"]::Vector{Float64}` (length K)
# When `per_dataset=true`, the metadata vectors are flat across datasets in
# the same dataset order as `_group_by_dataset` and are matched to clusters
# via `(dataset, id)` semantics — the n_clusters[d]-th entries belong to
# dataset d, etc. To keep this simple, the metadata is a single flat vector
# concatenated in dataset order.
#
# Algorithm (Campello 2013):
#   1. Core distances:  d_core(i) = distance to the min_pts-th NN of i.
#   2. Mutual reach.:   d_mreach(i,j) = max(d_core(i), d_core(j), euclidean(i,j))
#   3. MST on the complete d_mreach graph. Approximated by sparse k'-NN graph
#      with k' = knn_graph_k. Errors out if the resulting MST is disconnected.
#   4. Single-linkage hierarchy from MST.
#   5. Condense hierarchy: split events with both branches >= min_cluster_size
#      become parent->subcluster edges in the condensed tree; otherwise the
#      smaller branch's points "fall out" of the parent at that lambda.
#   6. Cluster persistence (stability): for each cluster C,
#        S(C) = Σ_{p falls out of C} (λ_fall(p, C) - λ_birth(C)) * size_at_fall
#   7. Cluster selection — :eom (default) or :leaf:
#        :eom  — walk bottom-up; keep cluster C if S(C) > Σ S(children selected),
#                else propagate child selections.
#        :leaf — select all leaves of the condensed tree.

"""
    HDBSCANConfig(; min_points=5, min_cluster_size=nothing, knn_graph_k=30,
                    cluster_selection_method=:eom, allow_single_cluster=false,
                    halo_trim_frac=0.10, use_3d=false, per_dataset=true,
                    remove_unclustered=false)

Configuration for HDBSCAN clustering of SMLM localizations. Pure-Julia
implementation of Campello/Moulavi/Sander 2013.

# Fields
- `min_points::Int = 5`: `k` for the core-distance computation. Larger values
  produce more conservative clusters.
- `min_cluster_size::Union{Int,Nothing} = nothing`: minimum cluster size in
  the condensed tree. When `nothing`, defaults to `min_points`.
- `knn_graph_k::Int = 30`: width `k'` of the sparse k'-NN graph used as the
  MST scaffold. Must be large enough that the resulting MST is connected;
  the backend errors out (with a clear message) if `knn_graph_k` is too small.
- `cluster_selection_method::Symbol = :eom`: how to pick clusters from the
  condensed tree — `:eom` (excess of mass; the canonical HDBSCAN choice) or
  `:leaf` (return all leaves of the condensed tree, like flat DBSCAN at
  varying ε).
- `allow_single_cluster::Bool = false`: when `true`, the root cluster (the
  whole connected mass) is a candidate for EOM selection. Useful for data
  that is essentially one tight blob with no internal real splits — the
  canonical EOM rule otherwise returns zero clusters in that case. Default
  matches the Python `hdbscan` package default.
- `halo_trim_frac::Float64 = 0.10`: halo-trim fraction. A point that fell out of its
  cluster individually (a transient / halo point) is kept only if it survived past
  this fraction of the cluster's λ-life; weakly-attached points that peeled off near
  the cluster's birth (density-connected background the MST routed onto the cluster's
  side) are dropped to noise. `0` disables trimming (raw condensed-tree labels — every
  point under a selected branch is a member); larger values trim more aggressively.
  Calibrate against the physical edge (the radius where cluster density crosses
  background); for the validated default the diffuse-cluster member radius matches
  that edge.
- `use_3d::Bool = false`: cluster in (x, y, z); requires 3D emitters.
- `per_dataset::Bool = true`: cluster within each `dataset` index independently.
- `remove_unclustered::Bool = false`: drop noise emitters from the output SMLD.

# Outputs
`cluster(smld, ::HDBSCANConfig)` returns `(smld_out, info)` matching the
shared `cluster` interface. Per-cluster *persistence* (stability) and birth
lambda are written to `smld_out.metadata` under the keys
`"hdbscan_cluster_persistence"` and `"hdbscan_cluster_lambda_birth"` —
flat vectors of length `info.n_clusters` in cluster-id order (and
concatenated across datasets when `per_dataset=true`).

# Example
```julia
cfg = HDBSCANConfig(min_points=10, min_cluster_size=20, knn_graph_k=50)
(smld_out, info) = cluster(smld, cfg)
persistence = smld_out.metadata["hdbscan_cluster_persistence"]
```

See also: [`AbstractClusterConfig`](@ref), [`ClusterInfo`](@ref), [`cluster`](@ref).
"""
Base.@kwdef struct HDBSCANConfig <: AbstractClusterConfig
    min_points::Int = 5
    min_cluster_size::Union{Int,Nothing} = nothing
    knn_graph_k::Int = 30
    cluster_selection_method::Symbol = :eom
    allow_single_cluster::Bool = false
    halo_trim_frac::Float64 = 0.10
    use_3d::Bool = false
    per_dataset::Bool = true
    remove_unclustered::Bool = false
end

# ----------------------------------------------------------------------------
# Internal: core HDBSCAN on a coordinate matrix (d × n, in microns by convention,
# but units pass through algebraically). Returns a NamedTuple with assignments,
# persistence, lambda_birth, n_clusters.
# ----------------------------------------------------------------------------
function _hdbscan_core(X::AbstractMatrix{Float64};
                      min_points::Int,
                      min_cluster_size::Int,
                      knn_graph_k::Int,
                      cluster_selection_method::Symbol,
                      allow_single_cluster::Bool = false,
                      halo_trim_frac::Float64 = 0.10)
    n = size(X, 2)
    if n == 0
        return (assignments = Int[],
                persistence = Float64[],
                lambda_birth = Float64[],
                n_clusters = 0)
    end
    if n < min_cluster_size
        # Not enough points to form any cluster.
        return (assignments = zeros(Int, n),
                persistence = Float64[],
                lambda_birth = Float64[],
                n_clusters = 0)
    end

    tree = NearestNeighbors.KDTree(X)

    # --- Step 1: core distances --------------------------------------------
    k_core = min(min_points, n - 1)
    core_dists = Vector{Float64}(undef, n)
    @inbounds for i in 1:n
        idxs, dists = NearestNeighbors.knn(tree, view(X, :, i), k_core + 1, true)
        # First entry is self at distance 0; take the (k_core+1)-th = k_core-th NN
        core_dists[i] = dists[end]
    end

    # --- Step 2/3: kNN graph with mutual reachability + Kruskal MST --------
    k_graph = min(knn_graph_k, n - 1)
    # Generous upper bound on edges; we dedup on the fly.
    edge_a = Int[]; edge_b = Int[]; edge_w = Float64[]
    sizehint!(edge_a, n * k_graph); sizehint!(edge_b, n * k_graph); sizehint!(edge_w, n * k_graph)
    @inbounds for i in 1:n
        idxs, dists = NearestNeighbors.knn(tree, view(X, :, i), k_graph + 1, true)
        for jj in 2:length(idxs)  # skip self
            j = idxs[jj]
            (i < j) || continue  # canonical orientation, dedup
            d = dists[jj]
            mreach = max(core_dists[i], core_dists[j], d)
            push!(edge_a, i); push!(edge_b, j); push!(edge_w, mreach)
        end
    end
    # Some near-pairs may only appear in one direction's kNN; add the reverse pairs we missed.
    @inbounds for i in 1:n
        idxs, dists = NearestNeighbors.knn(tree, view(X, :, i), k_graph + 1, true)
        for jj in 2:length(idxs)
            j = idxs[jj]
            (j < i) || continue
            # We need to make sure (j, i) is included — but the i-j edge from j's perspective
            # would already have been added in the previous loop. Skip; symmetry handled above.
        end
    end
    perm = sortperm(edge_w)

    # Union-find for Kruskal
    uf_parent = collect(1:n)
    uf_rank   = zeros(Int, n)
    function uf_find(x::Int)
        while uf_parent[x] != x
            uf_parent[x] = uf_parent[uf_parent[x]]
            x = uf_parent[x]
        end
        return x
    end

    mst_a = Int[]; mst_b = Int[]; mst_w = Float64[]
    sizehint!(mst_a, n - 1); sizehint!(mst_b, n - 1); sizehint!(mst_w, n - 1)
    @inbounds for k in perm
        a = edge_a[k]; b = edge_b[k]
        ra = uf_find(a); rb = uf_find(b)
        ra == rb && continue
        # union by rank
        if uf_rank[ra] < uf_rank[rb]
            ra, rb = rb, ra
        end
        uf_parent[rb] = ra
        if uf_rank[ra] == uf_rank[rb]; uf_rank[ra] += 1; end
        push!(mst_a, a); push!(mst_b, b); push!(mst_w, edge_w[k])
        length(mst_a) == n - 1 && break
    end
    # The kNN-graph approximation can leave the MST disconnected when
    # clusters are far apart. Repair by finding cross-component bridges:
    # for each point still in a multi-component state, query KDTree with
    # progressively larger k until we see a point in a different component.
    # Add the smallest such bridge per component, repeat until connected.
    while length(mst_a) < n - 1
        comp_id = [uf_find(i) for i in 1:n]
        bridges_a = Int[]; bridges_b = Int[]; bridges_w = Float64[]
        @inbounds for p in 1:n
            cp = comp_id[p]
            k_try = max(knn_graph_k, 2)
            found_b = 0; found_w = Inf
            while true
                k_use = min(k_try + 1, n)
                idxs, dists = NearestNeighbors.knn(tree, view(X, :, p), k_use, true)
                for jj in 2:length(idxs)
                    j = idxs[jj]
                    if comp_id[j] != cp
                        d = dists[jj]
                        mreach = max(core_dists[p], core_dists[j], d)
                        if mreach < found_w
                            found_b = j; found_w = mreach
                        end
                        @goto next_point
                    end
                end
                k_try == n - 1 && break
                k_try = min(k_try * 2, n - 1)
            end
            @label next_point
            if found_b > 0
                push!(bridges_a, p); push!(bridges_b, found_b); push!(bridges_w, found_w)
            end
        end
        isempty(bridges_a) && error("HDBSCAN: cannot bridge components — graph fundamentally disconnected.")
        bperm = sortperm(bridges_w)
        added_any = false
        @inbounds for k in bperm
            a = bridges_a[k]; b = bridges_b[k]
            ra = uf_find(a); rb = uf_find(b)
            ra == rb && continue
            if uf_rank[ra] < uf_rank[rb]; ra, rb = rb, ra; end
            uf_parent[rb] = ra
            if uf_rank[ra] == uf_rank[rb]; uf_rank[ra] += 1; end
            push!(mst_a, a); push!(mst_b, b); push!(mst_w, bridges_w[k])
            added_any = true
            length(mst_a) == n - 1 && break
        end
        added_any || error("HDBSCAN: bridge-step made no progress. Components: $(length(unique(comp_id))).")
    end

    # --- Step 4: single-linkage hierarchy ---------------------------------
    # Sort MST edges ascending (smallest weight first → highest lambda).
    mst_perm = sortperm(mst_w)
    n_internal = n - 1
    h_left   = Vector{Int}(undef, n_internal)
    h_right  = Vector{Int}(undef, n_internal)
    h_lambda = Vector{Float64}(undef, n_internal)
    h_size   = Vector{Int}(undef, n_internal)
    node_size = ones(Int, 2n - 1)  # leaves 1..n size=1; internal will be filled
    # Per-component union-find that tracks the current "top node id" for each component.
    uf2_parent = collect(1:n)
    uf2_top    = collect(1:n)        # current top tree node id (1..n at start)
    uf2_rank   = zeros(Int, n)
    function uf2_find(x::Int)
        while uf2_parent[x] != x
            uf2_parent[x] = uf2_parent[uf2_parent[x]]
            x = uf2_parent[x]
        end
        return x
    end
    @inbounds for (k, idx) in enumerate(mst_perm)
        a = mst_a[idx]; b = mst_b[idx]; w = mst_w[idx]
        ra = uf2_find(a); rb = uf2_find(b)
        left_id  = uf2_top[ra]
        right_id = uf2_top[rb]
        new_id = n + k
        sz = node_size[left_id] + node_size[right_id]
        node_size[new_id] = sz
        h_left[k]   = left_id
        h_right[k]  = right_id
        # lambda = 1/d. Guard against d == 0 (coincident points) → ∞.
        h_lambda[k] = w > 0 ? 1.0 / w : Inf
        h_size[k]   = sz
        # union by rank
        if uf2_rank[ra] < uf2_rank[rb]
            ra, rb = rb, ra
        end
        uf2_parent[rb] = ra
        uf2_top[ra]    = new_id
        if uf2_rank[ra] == uf2_rank[rb]; uf2_rank[ra] += 1; end
    end
    root_id = n + n_internal  # = 2n - 1

    # --- Step 5: condense hierarchy ---------------------------------------
    # Walk top-down from root. Maintain mapping from internal-node id to
    # condensed-cluster id, and emit (parent, child, lambda_fall, child_size)
    # records for each "fall out" event (both branch-becomes-subcluster and
    # branch-falls-as-noise cases — the difference is captured in child_size:
    # >= min_cluster_size means subcluster, == 1 means single-point fall,
    # in between means a chunk that fell out as a unit).
    #
    # Cluster id 1 = root (the trivial whole-data cluster).
    next_cluster_id = 1
    # Per-internal-node, the condensed cluster id it belongs to.
    cluster_of_internal = zeros(Int, n_internal)
    cluster_of_internal[n_internal] = 1  # root internal node belongs to root cluster
    cluster_birth_lambda = Float64[0.0]  # cluster 1 (root) birth = 0
    cluster_parent       = Int[0]        # cluster 1 has no parent
    cluster_children     = Vector{Vector{Int}}(); push!(cluster_children, Int[])

    # Leaf assignment + per-point fall-out lambda (lambda at which point left
    # its innermost cluster as an individual — only set if it fell out as
    # noise; otherwise the point gets the cluster's death lambda set later).
    point_cluster = ones(Int, n)
    point_lambda  = zeros(Float64, n)
    point_fell_individually = falses(n)

    # Condensed-tree records as 4 parallel vectors.
    ct_parent     = Int[]
    ct_child      = Int[]   # 0 = single point fall (multiple records, one per point), >0 = sub-cluster id
    ct_lambda     = Float64[]
    ct_child_size = Int[]

    # Iterative DFS using an explicit stack to avoid Julia's recursion limit at
    # 600k points. Each stack frame is a tuple (node_id, current_cluster).
    stack = Tuple{Int,Int}[]
    push!(stack, (root_id, 1))
    # We also need a separate "fall stack" for descend_falling traversals.
    # To keep one loop, encode mode: (node_id, cluster, lambda_falling) where
    # lambda_falling >= 0 means "falling" descent; -1.0 means "normal" descent.
    fall_stack = Tuple{Int,Int,Float64}[]

    # Drop a sub-tree of points into a cluster as either "transient" (they
    # left the cluster individually while it was still alive) or "stable"
    # (the cluster died and they were full members up to that moment).
    # Transient points are noise w.r.t. the cluster; stable ones inherit
    # the cluster's selected-ancestor label via walk-up.
    function descend_falling_iter!(start_node::Int, cluster::Int, lambda_fall::Float64;
                                  transient::Bool)
        empty!(fall_stack)
        push!(fall_stack, (start_node, cluster, lambda_fall))
        while !isempty(fall_stack)
            (nid, c, lam) = pop!(fall_stack)
            if nid <= n
                point_cluster[nid] = c
                point_lambda[nid]  = lam
                point_fell_individually[nid] = transient
                continue
            end
            k = nid - n
            push!(fall_stack, (h_left[k], c, lam))
            push!(fall_stack, (h_right[k], c, lam))
        end
    end

    # Walk top-down. Track each cluster's death lambda so we can set
    # point_lambda for points that stay until cluster death.
    cluster_death_lambda = Vector{Float64}()  # 1-indexed, parallel to other cluster_* vecs
    push!(cluster_death_lambda, 0.0)  # placeholder for root; set later if root ever splits

    while !isempty(stack)
        (nid, current_cluster) = pop!(stack)
        if nid <= n
            point_cluster[nid] = current_cluster
            # point_lambda left as 0 here; resolved post-walk via cluster death lambda
            continue
        end
        k = nid - n
        cluster_of_internal[k] = current_cluster
        L  = h_left[k]
        R  = h_right[k]
        lam = h_lambda[k]
        size_L = node_size[L]
        size_R = node_size[R]
        if size_L >= min_cluster_size && size_R >= min_cluster_size
            # Real split — current_cluster dies, two children born.
            cluster_death_lambda[current_cluster] = lam
            c_left  = next_cluster_id + 1
            c_right = next_cluster_id + 2
            next_cluster_id += 2
            push!(cluster_birth_lambda, lam)
            push!(cluster_birth_lambda, lam)
            push!(cluster_parent, current_cluster)
            push!(cluster_parent, current_cluster)
            push!(cluster_children, Int[])
            push!(cluster_children, Int[])
            push!(cluster_death_lambda, 0.0)
            push!(cluster_death_lambda, 0.0)
            push!(cluster_children[current_cluster], c_left)
            push!(cluster_children[current_cluster], c_right)
            # Condensed-tree records: both children fell out of parent at lam.
            push!(ct_parent, current_cluster); push!(ct_child, c_left);  push!(ct_lambda, lam); push!(ct_child_size, size_L)
            push!(ct_parent, current_cluster); push!(ct_child, c_right); push!(ct_lambda, lam); push!(ct_child_size, size_R)
            push!(stack, (L, c_left))
            push!(stack, (R, c_right))
        elseif size_L >= min_cluster_size
            # Right side is too small → transient fall out (cluster still alive).
            descend_falling_iter!(R, current_cluster, lam; transient = true)
            push!(ct_parent, current_cluster); push!(ct_child, 0); push!(ct_lambda, lam); push!(ct_child_size, size_R)
            push!(stack, (L, current_cluster))
        elseif size_R >= min_cluster_size
            descend_falling_iter!(L, current_cluster, lam; transient = true)
            push!(ct_parent, current_cluster); push!(ct_child, 0); push!(ct_lambda, lam); push!(ct_child_size, size_L)
            push!(stack, (R, current_cluster))
        else
            # Both sides too small — cluster dies. Points are stable members up to death.
            cluster_death_lambda[current_cluster] = lam
            descend_falling_iter!(L, current_cluster, lam; transient = false)
            descend_falling_iter!(R, current_cluster, lam; transient = false)
            push!(ct_parent, current_cluster); push!(ct_child, 0); push!(ct_lambda, lam); push!(ct_child_size, size_L)
            push!(ct_parent, current_cluster); push!(ct_child, 0); push!(ct_lambda, lam); push!(ct_child_size, size_R)
        end
    end

    # For points whose lambda wasn't set (they stayed in their innermost
    # cluster until that cluster died at split), set point_lambda to the
    # cluster's death lambda. If a cluster never died (e.g. a tiny dataset
    # where root is stable), we leave lambda at 0 — those points contribute
    # nothing to stability.
    @inbounds for p in 1:n
        if !point_fell_individually[p]
            c = point_cluster[p]
            point_lambda[p] = cluster_death_lambda[c]
        end
    end

    n_total_clusters = next_cluster_id
    # --- Step 6: stability per cluster ------------------------------------
    # Use the condensed tree records to credit each parent at the moment its
    # child fell out, weighted by child_size. This matches the canonical
    # HDBSCAN formula. Root's birth = 0 means root stability is unbounded
    # in principle; we exclude root from selection regardless (it's the
    # trivial whole-data cluster).
    stability = zeros(Float64, n_total_clusters)
    @inbounds for k in eachindex(ct_parent)
        p   = ct_parent[k]
        lam = ct_lambda[k]
        sz  = ct_child_size[k]
        b   = cluster_birth_lambda[p]
        if isfinite(lam) && isfinite(b)
            stability[p] += (lam - b) * sz
        end
    end

    # --- Step 7: cluster selection ----------------------------------------
    # Build child-list of clusters (already have cluster_children).
    # Selection arrays: selected[c] ∈ {true, false}; final cluster ids are 1..K
    # over the selected ones.
    selected = falses(n_total_clusters)
    if cluster_selection_method === :leaf
        # Select all leaves of the condensed tree (clusters with no condensed children).
        for c in 2:n_total_clusters  # skip root (cluster 1)
            isempty(cluster_children[c]) && (selected[c] = true)
        end
    elseif cluster_selection_method === :eom
        # Walk bottom-up. For each cluster, compare own stability to the
        # sum of its descendants' currently-best stabilities.
        # Order: process by decreasing cluster id (children get id > parent
        # by construction — we always birth two children with ids > parent).
        # Root (id 1) IS included as a selection candidate so that data with
        # one tight blob (no internal sub-splits) returns the whole mass as
        # one cluster instead of an empty result.
        best_stab = copy(stability)
        is_self_chosen = falses(n_total_clusters)
        # Bottom-up. The root (id=1) is included only when allow_single_cluster=true,
        # matching the Python hdbscan convention.
        bottom = allow_single_cluster ? 1 : 2
        for c in n_total_clusters:-1:bottom
            kids = cluster_children[c]
            if isempty(kids)
                best_stab[c] = stability[c]
                is_self_chosen[c] = true
            else
                child_sum = 0.0
                for ck in kids
                    child_sum += best_stab[ck]
                end
                if stability[c] >= child_sum
                    best_stab[c] = stability[c]
                    is_self_chosen[c] = true
                else
                    best_stab[c] = child_sum
                    is_self_chosen[c] = false
                end
            end
        end
        # Top-down selection.
        sel_stack = Int[]
        if allow_single_cluster
            push!(sel_stack, 1)
        else
            for ck in cluster_children[1]
                push!(sel_stack, ck)
            end
        end
        while !isempty(sel_stack)
            c = pop!(sel_stack)
            if is_self_chosen[c]
                selected[c] = true
            else
                for ck in cluster_children[c]
                    push!(sel_stack, ck)
                end
            end
        end
    else
        error("HDBSCANConfig.cluster_selection_method must be :eom or :leaf (got $(cluster_selection_method))")
    end

    # --- Final assignment: each point → its deepest selected ancestor cluster.
    # Walk up from point_cluster[p] via cluster_parent, taking the first
    # selected cluster encountered. If none selected, point is noise (id=0).
    selected_ids = findall(selected)  # original cluster ids 1..n_total_clusters
    id_map = Dict{Int,Int}()
    for (k, c) in enumerate(selected_ids)
        id_map[c] = k
    end

    # Fraction of a cluster's λ-life a transient (individually-fallen) point must
    # have survived to count as a member rather than weakly-attached halo. A point
    # that peeled off near the cluster's birth (≈ density-connected background that
    # the MST routed onto the cluster's side) is dropped to noise; one that survived
    # deep into the cluster's life is a genuine edge member and is kept. (GLOSH-style;
    # a blanket transient→noise rule over-prunes, since a Gaussian cluster sheds most
    # of its members transiently.)
    halo_min_life_frac = halo_trim_frac

    assignments = Vector{Int}(undef, n)
    @inbounds for p in 1:n
        # Walk up from the point's innermost cluster to find the deepest selected
        # ancestor. If none selected, the point is noise.
        c = point_cluster[p]
        while c > 0 && !selected[c]
            c = cluster_parent[c]
        end
        if !(c > 0 && selected[c])
            assignments[p] = 0
            continue
        end
        # Halo pruning: a transient point that left its cluster very close to the
        # cluster's birth is weakly attached → noise (see descend_falling_iter!).
        if point_fell_individually[p]
            cf    = point_cluster[p]
            birth = cluster_birth_lambda[cf]
            death = cluster_death_lambda[cf]
            if death > birth && (point_lambda[p] - birth) < halo_min_life_frac * (death - birth)
                assignments[p] = 0
                continue
            end
        end
        assignments[p] = id_map[c]
    end

    final_persistence = [stability[c] for c in selected_ids]
    final_birth       = [cluster_birth_lambda[c] for c in selected_ids]
    return (assignments = assignments,
            persistence = final_persistence,
            lambda_birth = final_birth,
            n_clusters = length(selected_ids))
end

# ----------------------------------------------------------------------------
# SMLD-facing dispatch
# ----------------------------------------------------------------------------
function cluster(smld::SMLMData.BasicSMLD, cfg::HDBSCANConfig)
    t0 = time_ns()
    smld = SMLMData.BasicSMLD(deepcopy(smld.emitters), smld.camera,
                              smld.n_frames, smld.n_datasets,
                              deepcopy(smld.metadata))
    n_in = length(smld.emitters)

    cfg.min_points >= 1 ||
        throw(ArgumentError("HDBSCANConfig.min_points must be ≥ 1 (got $(cfg.min_points))"))
    cfg.knn_graph_k >= 1 ||
        throw(ArgumentError("HDBSCANConfig.knn_graph_k must be ≥ 1 (got $(cfg.knn_graph_k))"))
    cfg.cluster_selection_method in (:eom, :leaf) ||
        throw(ArgumentError("HDBSCANConfig.cluster_selection_method must be :eom or :leaf (got $(cfg.cluster_selection_method))"))
    (0.0 <= cfg.halo_trim_frac < 1.0) ||
        throw(ArgumentError("HDBSCANConfig.halo_trim_frac must be in [0, 1) (got $(cfg.halo_trim_frac))"))
    mcs = cfg.min_cluster_size === nothing ? cfg.min_points : cfg.min_cluster_size
    mcs >= 2 ||
        throw(ArgumentError("HDBSCANConfig.min_cluster_size must be ≥ 2 (got $mcs)"))

    groups = _group_by_dataset(smld, cfg.per_dataset)
    cluster_sizes = Int[]
    n_clustered = 0
    persistence_all = Float64[]
    lambda_birth_all = Float64[]

    for idxs in groups
        isempty(idxs) && continue
        sub = view(smld.emitters, idxs)
        X = _coords_matrix(sub, cfg.use_3d)
        res = _hdbscan_core(X;
                            min_points = cfg.min_points,
                            min_cluster_size = mcs,
                            knn_graph_k = cfg.knn_graph_k,
                            cluster_selection_method = cfg.cluster_selection_method,
                            allow_single_cluster = cfg.allow_single_cluster,
                            halo_trim_frac = cfg.halo_trim_frac)
        @inbounds for (j, i) in pairs(idxs)
            smld.emitters[i].id = res.assignments[j]
        end
        # Cluster sizes for this group, indexed 1..n_clusters_local.
        if res.n_clusters > 0
            sizes_local = zeros(Int, res.n_clusters)
            @inbounds for a in res.assignments
                a > 0 && (sizes_local[a] += 1)
            end
            append!(cluster_sizes, sizes_local)
            n_clustered += sum(sizes_local)
            append!(persistence_all, res.persistence)
            append!(lambda_birth_all, res.lambda_birth)
        end
    end

    n_clusters = length(cluster_sizes)
    n_noise = n_in - n_clustered
    smld_out = _build_output(smld, cfg.remove_unclustered)
    # Stamp HDBSCAN-specific outputs onto metadata.
    smld_out.metadata["hdbscan_cluster_persistence"]  = persistence_all
    smld_out.metadata["hdbscan_cluster_lambda_birth"] = lambda_birth_all

    info = ClusterInfo(
        n_in,
        n_clustered,
        n_noise,
        n_clusters,
        cluster_sizes,
        :hdbscan,
        (time_ns() - t0) / 1e9,
    )
    return smld_out, info
end
