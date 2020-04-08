function bessel_downward_recurrence!(j, j′, x⁻¹::T, sinc, cosc, nmax, cf1, s;
                                     tol=100eps(T), verbose=false) where T
    big = √(typemax(T))
    invbig = inv(big)
    nj = length(j)
    if nmax > 1
        jₙ = s ? 1 : -1
        j′ₙ = cf1*jₙ

        S = (nmax-1)*x⁻¹
        for n = nmax:-1:2
            jₙ₋₁ = (S+x⁻¹)*jₙ + j′ₙ
            S -= x⁻¹
            j′ₙ₋₁ = S*jₙ₋₁ - jₙ

            jₙ = jₙ₋₁
            j′ₙ = j′ₙ₋₁

            if n-1 ≤ nj
                j[n-1] = jₙ₋₁
                j′[n-1] = j′ₙ₋₁
            end
        end
    end

    j′₀ = cosc - sinc*x⁻¹

    if nj > 1
        verbose && @show sinc
        ω = if abs(sinc) > 1e-1 # √(tol)
            sinc/j[1]
        else
            -j′₀/j[2]
        end
        lmul!(ω, j)
        lmul!(ω, j′)
    end

    j[1] = sinc
    j′[1] = j′₀
end

bessel_downward_recurrence!(::Nothing, args...; _...) = nothing

function neumann_upward_recurrence!(y, y′, x⁻¹::T, sinc, cosc; verbose=false, _...) where T
    y[1] = -cosc
    y′[1] = sinc + cosc*x⁻¹

    S = zero(T)
    for n = 2:length(y)
        y[n] = S*y[n-1] - y′[n-1]
        S += x⁻¹
        y′[n] = y[n-1] - (S+x⁻¹)*y[n]
    end
end

neumann_upward_recurrence!(::Nothing, args...; _...) = nothing
