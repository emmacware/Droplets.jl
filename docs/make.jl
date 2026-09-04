push!(LOAD_PATH, "../")
using Documenter, Droplets


Coalescence = [
        "SDM" => "coalescence.md",
        "Kernels" => "kernels.md",
]
Condensation = [
        "Condensation" => "condensation.md",
]

makedocs(
    sitename = "Droplets.jl",
    modules = [Droplets],
    format = Documenter.HTML(),
    pages = [
        "Home" => "index.md",
        "Coalescence" => Coalescence,
        "Condensation" => Condensation,
        "Turbulence" => "turbulence.md",
        "Advection" => "advection.md",
        "Thermodynamics" => "thermodynamics.md",
        "Single-Column Model" => "scm.md",
        "Droplet Attributes" => "droplets.md",
        "Settings" => "settings.md",
        "Constants" => "constants.md",
    ],
    clean = true,
)

deploydocs(repo = "github.com/emmacware/Droplets.jl.git",branch = "gh-pages", target = "build",forcepush=true)
