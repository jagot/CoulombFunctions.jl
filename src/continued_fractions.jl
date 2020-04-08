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
                        verbose=false) where T
    val_or_ϵ(x) = abs(x) < ϵ ? ϵ : x

    f = val_or_ϵ(b₀)
    C = f
    D = zero(T)
    δ = zero(T)

    fmt = if verbose
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
        if D < 0
            s = !s
        end
        Δ = C*D
        fn = f*Δ
        δ = abs(Δ-1)
        verbose && printfmtln(fmt, n, aₙ, bₙ, C, D, δ, fn, fn-f)
        f = fn
        δ < tol && return f,n,δ,s
    end

    @warn "Lentz–Thompson did not converge in $(max_iter) iterations, |Δ-1| = $(δ)"
    f,max_iter,δ,s
end

@doc raw"""
    bessel_fraction(x, n)

Evaluated the continued fraction for the spherical Bessel function

```math
\frac{j'_n(x)}{j_n(x)} =
\frac{n}{x} -
\frac{1}{T_{n+1}(x)-}\frac{1}{T_{n+2}(x)-}...\frac{1}{T_k-...},
```
where
```math
T_k(x) = \frac{2k+1}{x}.
```
"""
bessel_fraction(x::T, n::Integer; kwargs...) where T =
    lentz_thompson(n/x, k -> -one(T), k -> (2(n+k)+1)/x; kwargs...)
