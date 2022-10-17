using SphericalBesselFunctions
using Test

include("reference.jl")

@testset "SphericalBesselFunctions.jl" begin
    @testset "Accuracy" begin
        nx = 1001
        x = 10 .^ range(-1, stop=4, length=nx)
        nℓ = 105
        j, j′, y, y′ = @time bessels(x, nℓ)

        sel = 1:min(nℓ,104)

        nℓref = min(nℓ,length(sel))
        jref = zeros(nx, nℓref)
        j′ref = zeros(nx, nℓref)
        yref = zeros(nx, nℓref)
        y′ref = zeros(nx, nℓref)

        @time for (i,x) in enumerate(x)
            bessels_ref!(view(jref, i, :),
                         view(j′ref, i, :),
                         view(yref, i, :),
                         view(y′ref, i, :),
                         x)
        end

        j = j[:,sel]
        j′ = j′[:,sel]
        y = y[:,sel]
        y′ = y′[:,sel]
        jref = jref[:,sel]
        j′ref = j′ref[:,sel]
        yref = yref[:,sel]
        y′ref = y′ref[:,sel]

        @test j ≈ jref
        @test j′ ≈ j′ref

        # The irregular functions disagree wildly for x < ℓ, so we
        # mask that out
        mask = x .< (0:nℓ-2)'
        y[mask] .= 0
        y′[mask] .= 0
        yref[mask] .= 0
        y′ref[mask] .= 0

        # Just to make sure we did not mask out the whole matrix...
        @test any(!iszero, y)
        @test any(!iszero, y′)

        @test y ≈ yref
        @test y′ ≈ y′ref
    end

    @testset "Compare with tabulated data" begin
        j,j′,y,y′ = bessels(reference_data.z, reference_data.nℓ)
        @test j ≈ transpose(reference_data.j)
        @test j′ ≈ transpose(reference_data.j′)
    end
end
