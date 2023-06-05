using SmoothedParticleHydrodynamics
using Documenter

DocMeta.setdocmeta!(SmoothedParticleHydrodynamics, :DocTestSetup, :(using SmoothedParticleHydrodynamics); recursive=true)

makedocs(;
    modules=[SmoothedParticleHydrodynamics],
    authors="Alexander Barth <barth.alexander@gmail.com> and contributors",
    repo="https://github.com/Alexander-Barth/SmoothedParticleHydrodynamics.jl/blob/{commit}{path}#{line}",
    sitename="SmoothedParticleHydrodynamics.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://Alexander-Barth.github.io/SmoothedParticleHydrodynamics.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/Alexander-Barth/SmoothedParticleHydrodynamics.jl",
    devbranch="main",
)
