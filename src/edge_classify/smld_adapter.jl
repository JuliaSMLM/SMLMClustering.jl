"""
    classify_emitters(smld; kwargs...) -> EdgeClassificationResult

SMLD adapter: extracts `x_um`, `y_um` from `smld.emitters[].x/.y` and
`fov_um` from `smld.camera.pixel_edges_x/y`, then dispatches to the
coordinate-based core API.

If `smld_path` (a `String`) is provided instead of an SMLD object, the
file is loaded with `JLD2.load(smld_path, "smld")` and provenance metadata
(path, mtime, size) is recorded into `params.json`.
"""
function classify_emitters(smld; kwargs...)
    n = length(smld.emitters)
    x_um = Vector{Float64}(undef, n); y_um = Vector{Float64}(undef, n)
    @inbounds for i in 1:n
        x_um[i] = smld.emitters[i].x
        y_um[i] = smld.emitters[i].y
    end
    fov_um = (Float64(smld.camera.pixel_edges_x[1]),
              Float64(smld.camera.pixel_edges_x[end]),
              Float64(smld.camera.pixel_edges_y[1]),
              Float64(smld.camera.pixel_edges_y[end]))
    return classify_emitters(x_um, y_um;
                             fov_um = fov_um, kwargs...)
end
