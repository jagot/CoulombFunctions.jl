using Documenter, SphericalBesselFunctions

isdefined(Main, :NOPLOTS) && NOPLOTS || include("plots.jl")

makedocs(;
    modules=[SphericalBesselFunctions],
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://jagot.github.io/SphericalBesselFunctions.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Spherical Bessel Functions" => "spherical_bessel_functions.md",
        "Coulomb Functions" => "coulomb_functions.md",
        "Continued Fractions" => "continued_fractions.md",
    ],
    repo="https://github.com/jagot/SphericalBesselFunctions.jl/blob/{commit}{path}#L{line}",
    sitename="SphericalBesselFunctions.jl",
    authors="Stefanos Carlström <stefanos.carlstrom@gmail.com>",
    doctest=false,
)

deploydocs(;
    repo="github.com/jagot/SphericalBesselFunctions.jl",
    push_preview = true,
)
