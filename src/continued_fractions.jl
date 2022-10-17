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
                        ϵ = eps(real(T)), tol=100eps(real(T)),
                        verbosity=0) where T
    val_or_ϵ(x) = abs(x) < ϵ ? ϵ : x

    f = val_or_ϵ(b₀)
    C = f
    U = promote_type(typeof(a(1)), typeof(b(1)))
    D = zero(U)
    δ = zero(U)

    fmt = if verbosity>1
        if U <: Real
            printfmtln("{1:>4s} {2:>10s} {3:>10s} {4:>10s} {5:>10s} {6:>10s} {7:>10s} {8:>10s}",
                       "n", "a", "b", "C", "D", "|Δ-1|", "f", "δf")
            FormatExpr("{1:>4d} {2:10.3e} {3:10.3e} {4:10.3e} {5:10.3e} {6:10.3e} {7:10.3e} {8:10.3e}")
        else
            printfmtln("{1:>4s} {2:>10s} {3:>10s} {4:>10s} {5:>10s} {6:>10s} {7:>10s} {8:>10s} {9:>10s} {10:>10s} {11:>10s} {12:>10s} {13:>10s} {14:>10s}",
                       "n", "Re(a)", "Im(a)", "Re(b)", "Im(b)", "Re(C)", "Im(C)", "Re(D)", "Im(D)", "|Δ-1|", "Re(f)", "Im(f)", "Re(δf)", "Im(δf)")
            FormatExpr("{1:>4d} {2:10.3e} {3:10.3e} {4:10.3e} {5:10.3e} {6:10.3e} {7:10.3e} {8:10.3e} {9:10.3e} {10:10.3e} {11:10.3e} {12:10.3e} {13:10.3e} {14:10.3e}")
        end
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
        if verbosity>1
            if U <: Real
                printfmtln(fmt, n, aₙ, bₙ, C, D, δ, fn, fn-f)
            else
                printfmtln(fmt, n, real(aₙ), imag(aₙ), real(bₙ), imag(bₙ), real(C), imag(C), real(D), imag(D), δ, real(fn), imag(fn), real(fn-f), imag(fn-f))
            end
        end
        f = fn
        δ < 1e-6eps(abs(f)) && return f,n,δ,s,true
    end

    verbosity > 0 && @warn "Lentz–Thompson did not converge in $(max_iter) iterations, |Δ-1| = $(δ)"
    f,max_iter,δ,s,false
end

mutable struct Kahan{T}
    h::T
    hc::T
end

Kahan(h::T) where T = Kahan(h, zero(T))

function add!(k::Kahan, δh)
    y = δh - k.hc
    t = k.h + y
    δt = t - k.h
    k.hc = δt - y
    k.h = t
end

@doc raw"""
    steed_kahan(b₀, a, b)

Steed's algorithm for the forward evaluation of continued
fractions, with Kahan summation to avoid loss of accuracy.

```math
f_n = b_0 + \frac{a_1}{b_1+}\frac{a_2}{b_2+}...\frac{a_{n-1}}{b_{n-1}+}\frac{a_n}{b_n}
```

As described in

- Thompson, I., & Barnett, A. (1986). Coulomb and Bessel functions of
  complex arguments and order. Journal of Computational Physics,
  64(2), 490–509. http://dx.doi.org/10.1016/0021-9991(86)90046-x

"""
function steed_kahan(b₀::T, a, b;
                     max_iter=20000,
                     ϵ = eps(real(T)), tol=100eps(real(T)),
                     verbosity=0) where T
    val_or_ϵ(x) = abs(x) < ϵ ? ϵ : x

    D = inv(b(1))
    δh = a(1)*D

    U = promote_type(typeof(a(1)), typeof(b(1)))
    kahan = Kahan(U(b₀))
    h = add!(kahan, δh)

    fmt = if verbosity>1
        if U <: Real
            printfmtln("{1:>4s} {2:>10s} {3:>10s} {4:>10s} {5:>10s} {6:>10s}",
                       "n", "a", "b", "h", "D", "δh")
            FormatExpr("{1:>4d} {2:10.3e} {3:10.3e} {4:10.3e} {5:10.3e} {6:10.3e}")
        else
            # printfmtln("{1:>4s} {2:>10s} {3:>10s} {4:>10s} {5:>10s} {6:>10s} {7:>10s} {8:>10s} {9:>10s} {10:>10s} {11:>10s} {12:>10s} {13:>10s} {14:>10s}",
            #            "n", "Re(a)", "Im(a)", "Re(b)", "Im(b)", "Re(C)", "Im(C)", "Re(D)", "Im(D)", "|Δ-1|", "Re(f)", "Im(f)", "Re(δf)", "Im(δf)")
            # FormatExpr("{1:>4d} {2:10.3e} {3:10.3e} {4:10.3e} {5:10.3e} {6:10.3e} {7:10.3e} {8:10.3e} {9:10.3e} {10:10.3e} {11:10.3e} {12:10.3e} {13:10.3e} {14:10.3e}")
        end
    else
        nothing
    end

    s = true

    for n = 2:max_iter
        aₙ = a(n)
        bₙ = b(n)
        D = inv(val_or_ϵ(bₙ + aₙ*D))
        δh *= (bₙ*D-1)
        h = add!(kahan, δh)
        if isreal(D) && D < 0
            # See section 3.2 of
            #
            # - Barnett, A., Feng, D., Steed, J., & Goldfarb, L. (1974). Coulomb
            #   wave functions for all real $\eta$ and ρ. Computer Physics
            #   Communications, 8(5),
            #   377–395. http://dx.doi.org/10.1016/0010-4655(74)90013-7
            s = !s
        end
        if verbosity>1
            if U <: Real
                printfmtln(fmt, n, aₙ, bₙ, h, D, δh)
            # else
            #     printfmtln(fmt, n, real(aₙ), imag(aₙ), real(bₙ), imag(bₙ), real(C), imag(C), real(D), imag(D), δ, real(fn), imag(fn), real(fn-f), imag(fn-f))
            end
        end
        abs(δh) < 1e-6eps(abs(h)) && return h,n,δh,s,true
    end

    verbosity > 0 && @warn "Steed–Kahan did not converge in $(max_iter) iterations, δh = $(δh)"
    h,max_iter,δh,s,false
end
