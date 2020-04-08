"""
    powneg1(m)

Returns an integer power of negative unity, i.e. ``(-)^n``.
"""
powneg1(m::Integer) = iseven(m) ? 1 : -1

@doc raw"""
    reflect!(g, g′, offset=0)

Affect the reflection symmetries
```math
j_n(-z) = (-)^n j_n(z), \qquad
y_n(-z) = (-)^{n+1}y_n(z),
```
where the ``+1`` can be added using `offset`.
"""
function reflect!(g, g′, offset=0)
    n = 0:length(g)-1
    for (i,n) in enumerate(n)
        pn = powneg1(n+offset)
        g[i] *= pn
        g′[i]
    end
end

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
function bessels!(j::J, j′::J, y::Y, y′::Y, x::T; tol=100eps(T), kwargs...) where {J,Y,T<:Number}
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
        end
        if !isnothing(y)
            y .= -T(Inf)
            y′ .= T(Inf)
        end
        return
    end

    reflect = x < zero(T)
    reflect && (x = -x)

    cf1,_,_,s = bessel_fraction(x, ℓmax; kwargs...)

    x⁻¹ = inv(x)
    sinx,cosx = sincos(x)
    sinc,cosc = sinx*x⁻¹,cosx*x⁻¹

    bessel_downward_recurrence!(j, j′, x⁻¹, sinc, cosc, ℓmax+1, cf1, s; kwargs...)
    neumann_upward_recurrence!(y, y′, x⁻¹, sinc, cosc; kwargs...)

    reflect && (reflect!(j, j′); reflect!(y, y′, 1))
end

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
