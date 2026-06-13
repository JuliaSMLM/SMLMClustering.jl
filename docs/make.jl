using SMLMClustering
using Documenter

DocMeta.setdocmeta!(SMLMClustering, :DocTestSetup, :(using SMLMClustering); recursive=true)

makedocs(;
    modules=[SMLMClustering],
    authors="klidke@unm.edu",
    sitename="SMLMClustering.jl",
    format=Documenter.HTML(;
        canonical="https://JuliaSMLM.github.io/SMLMClustering.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Edge Classification" => "edge_classify.md",
    ],
)

deploydocs(;
    repo="github.com/JuliaSMLM/SMLMClustering.jl",
    devbranch="main",
)
