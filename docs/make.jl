using C13NV
using Documenter

DocMeta.setdocmeta!(C13NV, :DocTestSetup, :(using C13NV); recursive=true)

makedocs(;
    modules=[C13NV],
    authors="Michael Goerz <mail@michaelgoerz.net> and contributors",
    sitename="C13NV.jl",
    format=Documenter.HTML(;
        canonical="https://goerz.github.io/C13NV.jl",
        edit_link="master",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/goerz/C13NV.jl",
    devbranch="master",
)
