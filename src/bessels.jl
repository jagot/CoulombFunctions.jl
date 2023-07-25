# * Properties

"""
    powneg1(m)

Returns an integer power of negative unity, i.e. ``(-)^n``.
"""
powneg1(m::Integer) = iseven(m) ? 1 : -1

@doc raw"""
    reflect!(g, g′, offset=0)

Affect the reflection symmetries
```math
\tag{DLMF10.47.14}
\begin{aligned}
j_n(-z) &= (-)^n j_n(z), &
y_n(-z) &= (-)^{n+1}y_n(z), \\
j_n'(-z) &= (-)^{n-1} j_n'(z), &
y_n'(-z) &= (-)^{n} y_n'(z),
\end{aligned}
```
where the ``+1`` can be added using `offset`.
"""
function reflect!(g, g′, offset=0)
    n = 0:length(g)-1
    for (i,n) in enumerate(n)
        pn = powneg1(n+offset)
        g[i] *= pn
        g′[i] *= -pn
    end
end

# * Continued fractions

@doc raw"""
    bessel_fraction(x, n)

Evaluate the continued fraction for the spherical Bessel function

```math
\frac{j'_n(x)}{j_n(x)} =
\frac{n}{x} -
\frac{1}{T_{n+1}(x)-}\frac{1}{T_{n+2}(x)-}...\frac{1}{T_k(x)-...},
```
where
```math
T_k(x) = \frac{2k+1}{x}.
```
"""
bessel_fraction(x::T, n::Integer; cf_algorithm=lentz_thompson,
                max_iter=max(1000, 2ceil(Int, abs(x))), kwargs...) where T =
    cf_algorithm(n/x, k -> -one(T), k -> (2(n+k)+1)/x; max_iter=max_iter, kwargs...)

# * Recurrences

@doc raw"""
    bessel_downward_recurrence!(j, j′, x⁻¹, sinc, cosc, nmax, cf1, s; large)

Given the logarithmic derivative `cf1=j′[end]/j[end]` (computed using
[`bessel_fraction`](@ref)), and the sign `s`, fill in all lower orders
using the _downward recurrence_
```math
\begin{aligned}
g_{n-1} &= S_{n+1} g_n + g_n', &
g_{n-1}' &= S_{n-1} g_{n-1} - g_n, &
S_n &= \frac{n}{x},
\end{aligned}
```
and normalize the whole series using the [analytically known
expressions](https://dlmf.nist.gov/10.49)
```math
\tag{DLMF10.49.3}
\begin{aligned}
j_0(z) &= \frac{\sin z}{z}, &
j_0'(z) &= \frac{\cos z}{z} - \frac{\sin z}{z^2}.
\end{aligned}
```

If at any point during the recurrence, `j[i]>large`, renormalize
`j′[i:end] ./= j[i]` and `j[i:end] ./= j[i]` to avoid overflow, and
the continue the recurrence. This is very useful for small ``|x|`` and
large ``n``.

`nmax` allows starting the recurrence for a higher ``n`` than for
which there is space allocated in `j` and `j′` (this can help
convergence for the terms which are actually of interest).
"""
function bessel_downward_recurrence!(j, j′, x⁻¹::T, sinc, cosc, nmax, cf1, s;
                                     tol=100eps(real(T)), verbosity=0, large=∛(floatmax(real(T))), kwargs...) where T
    verbosity > 0 && @show x⁻¹, sinc, cosc, nmax, cf1, s
    nj = length(j)
    nj == 0 && return
    if nmax > 1
        jₙ = s ? 1 : -1
        j′ₙ = cf1*jₙ

        S = (nmax-1)*x⁻¹
        @inbounds for n = nmax:-1:2
            jₙ₋₁ = (S+x⁻¹)*jₙ + j′ₙ
            S -= x⁻¹
            j′ₙ₋₁ = S*jₙ₋₁ - jₙ

            if abs(jₙ₋₁) > large
                verbosity > 1 && @info "Bessel downward recurrence larger than $(large), renormalizing" n jₙ₋₁ j′ₙ₋₁
                ij = inv(jₙ₋₁)
                jₙ₋₁ = one(jₙ₋₁)
                j′ₙ₋₁ *= ij
                lmul!(ij, view(j, n:nj))
                lmul!(ij, view(j′, n:nj))
            end

            jₙ = jₙ₋₁
            j′ₙ = j′ₙ₋₁

            if n-1 ≤ nj
                j[n-1] = jₙ₋₁
                j′[n-1] = j′ₙ₋₁
            end
        end
    end

    j′₀ = cosc - sinc*x⁻¹

    @inbounds if nj > 1
        verbosity > 0 && @show sinc
        ω = if abs(sinc) > 1e-1 # √(tol)
            sinc/j[1]
        else
            -j′₀/j[2]
        end
        verbosity > 0 && @show ω
        lmul!(ω, j)
        lmul!(ω, j′)
    end

    j[1] = sinc
    j′[1] = j′₀
end

bessel_downward_recurrence!(::Nothing, args...; _...) = nothing

@doc raw"""
    neumann_upward_recurrence(y, y′, x⁻¹, sinc, cosc)

Generate the irregular Neumann functions using (stable) upward
recurrence
```math
\begin{aligned}
g_{n+1} &= S_n g_n - g_n', &
g_{n+1}' &= g_n - S_{n+2} g_{n+1},
\end{aligned}
```
starting from the [analytically known
expressions](https://dlmf.nist.gov/10.49)
```math
\tag{DLMF10.49.5}
\begin{aligned}
y_0(z) &= -\frac{\cos z}{z}, &
y_0'(z) &= \frac{\sin z}{z} + \frac{\cos z}{z^2}.
\end{aligned}
```
"""
function neumann_upward_recurrence!(y, y′, x⁻¹::T, sinc, cosc; _...) where T
    y[1] = -cosc
    y′[1] = sinc + cosc*x⁻¹

    S = zero(T)
    @inbounds for n = 2:length(y)
        y[n] = S*y[n-1] - y′[n-1]
        S += x⁻¹
        y′[n] = y[n-1] - (S+x⁻¹)*y[n]
    end
end

neumann_upward_recurrence!(::Nothing, args...; _...) = nothing

# * Driver

@doc raw"""
    bessels!(j, j′, y, y′, x)

Compute all spherical Bessel (regular: ``j_n``) and Neumann
(irregular: ``y_n``) functions and their derivatives at `x` and store
the results in the vectors `j`, `j′`, `y`, `y′`. If only the Bessel
functions or the Neumann functions are of interest, the other pair of
arrays can be substituted by `nothing`. However, it is not possible to
compute only the functions but not the derivatives, since they are
generated using the following recurrence relations:

```math
\begin{aligned}
g_{n-1} &= \frac{n+1}{x} g_n + g_n'\\
g'_{n-1} &= \frac{n-1}{x} g_{n-1} - g_n.
\end{aligned}
```

These recurrence relations are employed in a _downward_ fashion for
the Bessel functions and an _upward_ fashion for the Neumann
functions.

It is assumed that all passed arrays are of the same lengths (not
checked).

"""
function bessels!(j::J, j′::J, y::Y, y′::Y, x::T; tol=100eps(real(T)), verbosity=0, kwargs...) where {J,Y,T<:Number}
    ℓmax = if isnothing(j)
        # No output requested
        isnothing(y) && return
        length(y)-1
    else
        nj = length(j)
        !isnothing(y) && !(length(y) == nj == length(j′) == length(y′)) &&
            throw(DimensionMismatch("Lengths of j and y and their derivatives must agree"))
        nj-1
    end

    if iszero(x)
        if !isnothing(j)
            j .= zero(T)
            j[1] = one(T)
            j′ .= zero(T)
            if ℓmax > 0
                j′[2] = one(T)/3 # Derivative of Eq. 10.52.1 https://dlmf.nist.gov/10.52
            end
        end
        if !isnothing(y)
            # This assumes that the limit is approached from positive x.
            y .= -T(Inf)
            y′ .= T(Inf)
        end
        return
    end

    reflect = real(x) < zero(real(T))
    reflect && (x = -x)

    cf1,_,_,s,converged = bessel_fraction(x, ℓmax+1; verbosity=verbosity-1, kwargs...)
    converged || verbosity > 0 && @info "Consider increasing ℓmax beyond $(ℓmax)"

    x⁻¹ = inv(x)
    sinx,cosx = sincos(x)
    sinc,cosc = sinx*x⁻¹,cosx*x⁻¹

    bessel_downward_recurrence!(j, j′, x⁻¹, sinc, cosc, ℓmax+2, cf1, s; verbosity=verbosity-2, kwargs...)
    neumann_upward_recurrence!(y, y′, x⁻¹, sinc, cosc; kwargs...)

    reflect && (reflect!(j, j′); reflect!(y, y′, 1))
end

# * Interface

"""
    bessels!(j, j′, y, y′, x::AbstractVector; kwargs...)

Loop through all values of `x` and compute all Bessel and Neumann
functions, storing the results in the preallocated matrices `j`, `j′`,
`y`, `y′`.
"""
function bessels!(j, j′, y, y′, x::AbstractVector; kwargs...)
    size(j,1) == size(j′,1) == size(y,1) == size(y′,1) == length(x) &&
        size(j,2) == size(j′,2) == size(y,2) == size(y′,2) ||
        throw(DimensionError("The dimension of the output arrays do not agree"))
    for (i,x) in enumerate(x)
        bessels!(view(j, i, :),
                 view(j′, i, :),
                 view(y, i, :),
                 view(y′, i, :),
                 x; kwargs...)
    end
end

"""
    bessels(x::AbstractVector, nℓ; kwargs...)

Convenience wrapper around [`bessels!`](@ref) that preallocates output
matrices of the appropriate dimensions.
"""
function bessels(x::AbstractVector{T}, nℓ; kwargs...) where T
    nx = length(x)
    j = zeros(T, nx, nℓ)
    j′ = zeros(T, nx, nℓ)
    y = zeros(T, nx, nℓ)
    y′ = zeros(T, nx, nℓ)
    bessels!(j, j′, y, y′, x; kwargs...)
    j, j′, y, y′
end

export bessels!, bessels
