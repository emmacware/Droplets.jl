using Documenter, Droplets

makedocs(
    sitename = "Docs for Droplets.jl",
    doctest = false,
    format = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true"),
    remotes = nothing,
    clean = false,
    pages = Any[
        "Home" => "index.md",
    ],
)

deploydocs(repo = "github.com/emmacware/droplets.jl", target = "build")
