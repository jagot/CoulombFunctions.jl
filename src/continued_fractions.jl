@doc raw"""
    lentz_thompson(b₀, a, b)

Lentz–Thompson algorithm for the forward evaluation of continued
fractions:

```math
f_n = b_0 + \frac{a_1}{b_1+}\frac{a_2}{b_2+}...\frac{a_{n-1}}{b_{n-1}+}\frac{a_n}{b_n}
```

As described in

- Barnett, A. R. (1996). The Calculation of Spherical Bessel and
  Coulomb Functions. In (Eds.), Computational Atomic Physics
  (pp. 181–202). Springer Berlin Heidelberg.

"""
function lentz_thompson(b₀::T, a, b;
                        max_iter=20000,
                        ϵ = eps(T), tol=100eps(T),
                        verbosity=0) where T
    val_or_ϵ(x) = abs(x) < ϵ ? ϵ : x

    f = val_or_ϵ(b₀)
    C = f
    D = zero(T)
    δ = zero(T)

    fmt = if verbosity>1
        printfmtln("{1:>4s} {2:>10s} {3:>10s} {4:>10s} {5:>10s} {6:>10s} {7:>10s} {8:>10s}",
                   "n", "a", "b", "C", "D", "|Δ-1|", "f", "δf")
        FormatExpr("{1:>4d} {2:10.3e} {3:10.3e} {4:10.3e} {5:10.3e} {6:10.3e} {7:10.3e} {8:10.3e}")
    else
        nothing
    end

    s = true

    for n = 1:max_iter
        aₙ = a(n)
        bₙ = b(n)
        C = val_or_ϵ(bₙ + aₙ/C)
        D = inv(val_or_ϵ(bₙ + aₙ*D))
        if isreal(D) && D < 0
            # See section 3.2 of
            #
            # - Barnett, A., Feng, D., Steed, J., & Goldfarb, L. (1974). Coulomb
            #   wave functions for all real $\eta$ and ρ. Computer Physics
            #   Communications, 8(5),
            #   377–395. http://dx.doi.org/10.1016/0010-4655(74)90013-7
            s = !s
        end
        Δ = C*D
        fn = f*Δ
        δ = abs(Δ-1)
        verbosity>1 && printfmtln(fmt, n, aₙ, bₙ, C, D, δ, fn, fn-f)
        f = fn
        δ < tol && return f,n,δ,s,true
    end

    verbosity > 0 && @warn "Lentz–Thompson did not converge in $(max_iter) iterations, |Δ-1| = $(δ)"
    f,max_iter,δ,s,false
end
