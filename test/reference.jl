using SpecialFunctions
using LinearAlgebra

using DelimitedFiles

function bessels_ref!(j, j′, y, y′, x)
    for i = 1:length(j)
        n = i - 1
        ν = (2n + 1)/2
        j[i] = besselj(ν, x)
        y[i] = bessely(ν, x)
    end

    ω = √(π/2x)
    lmul!(ω, j)
    lmul!(ω, y)

    j′[1] = -y[1]-j[1]/x
    y′[1] = j[1]-y[1]/x

    for i = 2:length(j)
        n = i - 1
        ν = (n + 1)/x
        j′[i] = j[i-1] - ν*j[i]
        y′[i] = y[i-1] - ν*y[i]
    end
end

function load_bessel_ref(filename)
    d = readdlm(filename)
    nℓ = (size(d,1)-1)÷2

    z = vec(d[1,:])
    j = d[1 .+ (1:nℓ),:]
    j′ = d[1 + nℓ .+ (1:nℓ),:]

    (nℓ=nℓ, z=z, j=j, j′=j′)
end

bessel_reference_data_1 = load_bessel_ref(joinpath(dirname(@__FILE__), "bessel-ref.txt"))
bessel_reference_data_2 = load_bessel_ref(joinpath(dirname(@__FILE__), "bessel-ref-L=1000.txt"))
