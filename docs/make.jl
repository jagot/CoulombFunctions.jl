using Documenter, SphericalBesselFunctions

isdefined(Main, :NOPLOTS) && NOPLOTS || include("plots.jl")

makedocs(;
    modules=[SphericalBesselFunctions],
    format=Documenter.HTML(),
    pages=[
        "Home" => "index.md",
    ],
    repo="https://github.com/jagot/SphericalBesselFunctions.jl/blob/{commit}{path}#L{line}",
    sitename="SphericalBesselFunctions.jl",
    authors="Stefanos Carlström <stefanos.carlstrom@gmail.com>"
)

deploydocs(;
    repo="github.com/jagot/SphericalBesselFunctions.jl",
)
