using C13NV
using Documenter
import Pkg


PROJECT_TOML = Pkg.TOML.parsefile(joinpath(@__DIR__, "..", "Project.toml"))
VERSION = PROJECT_TOML["version"]
NAME = PROJECT_TOML["name"]
AUTHORS = join(PROJECT_TOML["authors"], ", ") * " and contributors"
GITHUB = "https://github.com/goerz/C13NV.jl"

println("Starting makedocs")

PAGES = ["Home" => "index.md",]

makedocs(;
    authors = "Michael Goerz <mail@michaelgoerz.net> and contributors",
    sitename = "C13NV.jl",
    format = Documenter.HTML(;
        prettyurls = true,
        canonical = "https://goerz.github.io/C13NV.jl",
        edit_link = "master",
        footer = "[$NAME.jl]($GITHUB) v$VERSION docs powered by [Documenter.jl](https://github.com/JuliaDocs/Documenter.jl).",
        assets = String[],
    ),
    pages = PAGES,
)

println("Finished makedocs")

deploydocs(; repo = "github.com/goerz/C13NV.jl", devbranch = "master", push_preview = true)
