using SMLMClustering
using Documenter

DocMeta.setdocmeta!(SMLMClustering, :DocTestSetup, :(using SMLMClustering); recursive=true)

makedocs(;
    modules = [SMLMClustering, SMLMClustering.EdgeClassify],
    authors = "klidke@unm.edu",
    sitename = "SMLMClustering.jl",
    format = Documenter.HTML(;
        # github.io hosts are lower-cased
        canonical = "https://juliasmlm.github.io/SMLMClustering.jl",
        edit_link = "main",
        assets = String[],
    ),
    pages = [
        "Home" => "index.md",
        "User Guide" => "guide.md",
        "Methods" => [
            "Overview" => "methods/index.md",
            "DBSCAN" => "methods/dbscan.md",
            "HDBSCAN" => "methods/hdbscan.md",
            "Hierarchical" => "methods/hierarchical.md",
            "Voronoi (SR-Tesseler)" => "methods/voronoi.md",
            "MRF density-regime" => "methods/mrf_density.md",
            "Point hysteresis" => "methods/point_hysteresis.md",
            "Hopkins statistic" => "methods/hopkins.md",
            "Voronoi density" => "methods/voronoi_density.md",
            "Local contrast" => "methods/local_contrast.md",
            "Edge classification" => "methods/edge_classify.md",
        ],
        "Results Gallery" => "results.md",
        "API Reference" => "api.md",
    ],
    # checkdocs=:exports enforces every exported symbol is documented; warnonly is
    # temporary during the docs build-out so a missing docstring / cross-ref surfaces
    # as a warning instead of failing the build. Tighten to `warnonly = false` before release.
    checkdocs = :exports,
    warnonly = true,
)

deploydocs(;
    repo = "github.com/JuliaSMLM/SMLMClustering.jl",
    devbranch = "main",
)
