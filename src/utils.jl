# Shared internal helpers — coordinate extraction and distance computation.
# Loaded before backend files; all functions are package-private (underscore prefix).

using Clustering

# Build a d×n coordinate matrix in microns from a vector of emitters.
# d = 2 when use_3d=false; d = 3 when use_3d=true (requires :z property).
function _coords_matrix(emitters::AbstractVector{<:SMLMData.AbstractEmitter}, use_3d::Bool)
    n = length(emitters)
    if use_3d
        isempty(emitters) || hasproperty(first(emitters), :z) ||
            error("use_3d=true requires 3D emitters (e.g. Emitter3DFit); got $(eltype(emitters)).")
        X = Matrix{Float64}(undef, 3, n)
        @inbounds for i in 1:n
            e = emitters[i]
            X[1, i] = e.x
            X[2, i] = e.y
            X[3, i] = e.z
        end
        return X
    else
        X = Matrix{Float64}(undef, 2, n)
        @inbounds for i in 1:n
            e = emitters[i]
            X[1, i] = e.x
            X[2, i] = e.y
        end
        return X
    end
end

# Symmetric n×n pairwise Euclidean distance matrix from a d×n column-major matrix.
function _pairwise_distances(X::Matrix{Float64})
    d, n = size(X)
    D = Matrix{Float64}(undef, n, n)
    @inbounds for j in 1:n
        D[j, j] = 0.0
        for i in (j + 1):n
            dist = 0.0
            for k in 1:d
                diff = X[k, i] - X[k, j]
                dist += diff * diff
            end
            dist = sqrt(dist)
            D[i, j] = dist
            D[j, i] = dist
        end
    end
    D
end
