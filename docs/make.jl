# run me using `julia make.jl` in this directory.

using Documenter, Clusters
makedocs(
sitename="Clusters.jl",
modules=[Clusters],
format=Documenter.HTML(prettyurls = get(ENV, "CI", nothing)=="true"),
pages=["Home" => "index.md"],
strict=true,
)

deploydocs(repo="github.com/michael-petersen/Clusters.jl",
            push_preview=true)
