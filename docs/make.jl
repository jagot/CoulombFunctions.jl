using Documenter, SphericalBesselFunctions

makedocs(;
    modules=[SphericalBesselFunctions],
    format=Documenter.HTML(),
    pages=[
        "Home" => "index.md",
    ],
    repo="https://github.com/jagot/SphericalBesselFunctions.jl/blob/{commit}{path}#L{line}",
    sitename="SphericalBesselFunctions.jl",
    authors="Stefanos Carlström <stefanos.carlstrom@gmail.com>",
    assets=String[],
)

deploydocs(;
    repo="github.com/jagot/SphericalBesselFunctions.jl",
)
