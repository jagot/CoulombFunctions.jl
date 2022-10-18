using SphericalBesselFunctions
import SphericalBesselFunctions: coulomb_fraction1, coulomb_fraction2, lentz_thompson, steed_kahan
using Test

include("reference.jl")

@testset "SphericalBesselFunctions.jl" begin
    @testset "Spherical Bessel functions" begin
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

        @testset "Compare regular Bessel functions, with tabulated data, ℓₘₐₓ=$(rd.nℓ-1)" for rd in [bessel_reference_data_1, bessel_reference_data_2]
            j,j′,y,y′ = bessels(rd.z, rd.nℓ)
            @test j ≈ transpose(rd.j)
            @test j′ ≈ transpose(rd.j′)
        end
    end

    @testset "Coulomb functions" begin
        @testset "Continued fractions" begin
            x = 10.0 .^ (-2:4)
            η = 10.0 .^ (-1:3)
            η = sort(vcat(-η,0,η))
            λ = [0,0.1,0.5,1,13.14,100,1000.6]

            nx = length(x)
            nη = length(η)
            nλ = length(λ)

            cf1_ref = load_cf_ref("coulomb-cf1-ref", nx, nη, nλ)
            cf2_ref = load_cf_ref("coulomb-cf2-real-ref", nx, nη, nλ) +
                im*load_cf_ref("coulomb-cf2-imag-ref", nx, nη, nλ)

            @testset "Coulomb continued fraction 1, $(label)" for (alg, label) in [(lentz_thompson, "Lentz–Thompson"),
                                                                                   (steed_kahan, "Steed–Kahan")]
                tol = √(eps())
                for (i,x) in enumerate(x)
                    for (j,η) in enumerate(η)
                        for (k,λ) in enumerate(λ)
                            cf1,iters,δf,s,isconverged = coulomb_fraction1(x, η, λ, cf_algorithm=alg)
                            @test cf1 ≈ cf1_ref[i,j,k]
                            @test δf < tol
                            @test isconverged
                        end
                    end
                end

                @testset "Verbose" begin
                    cf1,iters,δf,s,isconverged = coulomb_fraction1(1.0, -1.0, 100, cf_algorithm=alg, verbosity=4)
                    @test cf1 ≈ 100.98517182192248
                    @test δf < tol
                    @test isconverged
                end
            end
            @testset "Coulomb continued fraction 2, $(label)" for (alg, label) in [(lentz_thompson, "Lentz–Thompson"),
                                                                                   (steed_kahan, "Steed–Kahan")]
                tol = √(eps())
                for (i,x) in enumerate(x)
                    for (j,η) in enumerate(η)
                        for (k,λ) in enumerate(λ)
                            cf2,iters,δf,s,isconverged = coulomb_fraction2(x, η, λ, 1, cf_algorithm=alg)
                            @test cf2 ≈ cf2_ref[i,j,k]
                            @test δf < tol
                            @test isconverged
                        end
                    end
                end

                @testset "Verbose" begin
                    cf2,iters,δf,s,isconverged = coulomb_fraction2(1.0, -1.0, 100, 1, cf_algorithm=alg, verbosity=4)
                    @test cf2 ≈ -99.98497373591508 - 3.0677993260103584e-13im
                    @test δf < tol
                    @test isconverged
                end
            end
        end
    end
end
