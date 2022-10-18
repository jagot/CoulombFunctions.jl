using SpecialFunctions
using LinearAlgebra

using DelimitedFiles

using Test

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

ref_file(name) = joinpath(dirname(@__FILE__), name*".txt")
load_ref(name) = readdlm(ref_file(name))

function load_bessel_ref(name)
    d = load_ref(name)
    nℓ = (size(d,1)-1)÷2

    z = vec(d[1,:])
    j = d[1 .+ (1:nℓ),:]
    j′ = d[1 + nℓ .+ (1:nℓ),:]

    (nℓ=nℓ, z=z, j=j, j′=j′)
end

bessel_reference_data_1 = load_bessel_ref("bessel-ref")
bessel_reference_data_2 = load_bessel_ref("bessel-ref-L=1000")

load_cf_ref(name, sizes...) = reshape(load_ref(name), sizes...)

function compare_with_coulomb_reference(name, x, η, ℓs)
    @testset "Coulomb accuracy $(name)" begin
        f_ref = load_ref("coulomb-ref-"*name*"-f")
        f′_ref = load_ref("coulomb-ref-"*name*"-fprime")
        if isfile("coulomb-ref-"*name*"-g")
            g_ref = load_ref("coulomb-ref-"*name*"-g")
            g′_ref = load_ref("coulomb-ref-"*name*"-gprime")
            f,f′,g,g′ = coulombs(x, η, ℓs)

            @test f ≈ f_ref
            @test f′ ≈ f′_ref
            @test g ≈ g_ref
            @test g′ ≈ g′_ref
        else
            f,f′ = coulombFs(x, η, ℓs)

            @test f ≈ f_ref
            @test f′ ≈ f′_ref
        end
    end
end
