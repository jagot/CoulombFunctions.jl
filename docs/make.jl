using Documenter, CoulombFunctions

isdefined(Main, :NOPLOTS) && NOPLOTS || include("plots.jl")

makedocs(;
    modules=[CoulombFunctions],
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://jagot.github.io/CoulombFunctions.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Spherical Bessel Functions" => "spherical_bessel_functions.md",
        "Coulomb Functions" => "coulomb_functions.md",
        "Continued Fractions" => "continued_fractions.md",
    ],
    repo="https://github.com/jagot/CoulombFunctions.jl/blob/{commit}{path}#L{line}",
    sitename="CoulombFunctions.jl",
    authors="Stefanos Carlström <stefanos.carlstrom@gmail.com>",
    doctest=false,
)

deploydocs(;
    repo="github.com/jagot/CoulombFunctions.jl",
    push_preview = true,
)
