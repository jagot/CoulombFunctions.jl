import SpecialFunctions
const Γ = SpecialFunctions.gamma

# * Properties
@doc raw"""
    turning_point(η, ℓ)

Computes the turning point ``\rho_{\mathrm{TP}}`` of the Coulomb
functions, which separates the monotonic region
``\rho<\rho_{\mathrm{TP}}`` from the oscillatory region
``\rho>\rho_{\mathrm{TP}}``.

See Eq. (2) of

- Barnett, A. (1982). COULFG: Coulomb and Bessel Functions and Their
  Derivatives, for Real Arguments, by Steed's Method. Computer Physics
  Communications, 27(2),
  147–166. http://dx.doi.org/10.1016/0010-4655(82)90070-4

or

- Barnett, A. R. (1996). The Calculation of Spherical Bessel and
  Coulomb Functions. In (Eds.), Computational Atomic Physics
  (pp. 181–202). : Springer Berlin Heidelberg.
"""
turning_point(η, ℓ) = η + √(η^2 + ℓ*(ℓ+1))

# https://dlmf.nist.gov/33.2
coulomb_normalization(η, ℓ) = 2^ℓ*exp(-π*η/2)*abs(Γ(ℓ+1+im*η))/Γ(2ℓ+2)

# * Continued fractions

cf1R²(k,η) = 1 + η^2/k^2
cf1S(k,η,x) = k/x + η/k
cf1T(k,η,x) = (2k+1)*(inv(x) + η/(k^2 + k))

@doc raw"""
    coulomb_fraction1(x, η, n)

Evaluate the continued fraction for the logarithmic derivative of the
Coulomb function

```math
\frac{F'_n(\eta,x)}{F_n(\eta,x)} =
S_{n+1}(\eta,x) -
\frac{R_{n+1}^2(\eta)}{T_{n+1}(\eta,x)-}\frac{R_{n+2}^2(\eta)}{T_{n+2}(\eta,x)-}...\frac{R_k^2(\eta)}{T_k(\eta,x)-...},
```
where
```math
\begin{aligned}
R_k(\eta) &= \sqrt{1+\frac{\eta^2}{k^2}}, \qquad
S_k(\eta,x) = \frac{k}{x} + \frac{\eta}{k}, \\
T_k(\eta,x) &= S_k(\eta,x) + S_{k+1}(\eta,x) = (2k+1)\left(\frac{1}{x} + \frac{\eta}{k^2+k}\right).
\end{aligned}
```
"""
coulomb_fraction1(x::T, η::T, n::Integer; cf_algorithm=lentz_thompson,
                  max_iter=1000 + max(1, 5ceil(Int, √(abs(x^2-2η*x)))), kwargs...) where T =
    cf_algorithm(cf1S(n+1,η,x),
                 k -> -cf1R²(n+k,η),
                 k -> cf1T(n+k,η,x); max_iter=max_iter, kwargs...)

@doc raw"""
    coulomb_fraction2(x, η, n, ω)

Evaluate the continued fraction for the logarithmic derivative of the
Hankel function ``H^\omega_\lambda(\eta,x) = G_\lambda(\eta,x) +
\mathrm{i}\omega F_\lambda(\eta,x)``

```math
\frac{H^{\omega'}_n(\eta,x)}{H^\omega_n(\eta,x)} =
p+\mathrm{i}q =
\mathrm{i}\omega\left(1-\frac{\eta}{x}\right) +
\frac{\mathrm{i}\omega}{x}
\frac{ac}{2(x-\eta+\mathrm{i}\omega)+}
\frac{(a+1)(c+1)}{2(x-\eta+2\mathrm{i}\omega)+...},
```
where
```math
\begin{aligned}
a &= 1+n+\mathrm{i}\omega\eta, &
b &= 2n + 2, &
c &= -n + \mathrm{i}\omega\eta.
\end{aligned}
```
"""
function coulomb_fraction2(x::T, η::T, n::Integer, ω; cf_algorithm=lentz_thompson,
                           max_iter=max(1000, 2ceil(Int, 5000/abs(x))), kwargs...) where T
    imω = im*ω
    r = cf_algorithm(x-η,
                     k -> (imω*η - n - 1 + k)*(imω*η + n + k),
                     k -> 2*(x - η + imω*k); max_iter=max_iter, kwargs...)
    ((imω/x)*r[1],r[2:end]...)
end

# * Recurrences

function coulomb_downward_recurrence!(F, F′, x⁻¹::T, η, ℓ::UnitRange, cf1, s; verbosity=0,
                                      large=∛(floatmax(real(T))), kwargs...) where T
    Fₙ = s ? 1 : -1
    F′ₙ = cf1*Fₙ
    η² = η^2
    nF = length(F)
    for i = nF:-1:1
        n = ℓ[i]
        S = n*x⁻¹ + η/n
        R = √(1 + η²/n^2)

        Fₙ₋₁ = (S*Fₙ + F′ₙ)/R
        F′ₙ₋₁ = S*Fₙ₋₁ - R*Fₙ

        if abs(Fₙ₋₁) > large
            verbosity > 1 && @info "Coulomb downward recurrence larger than $(large), renormalizing" n Fₙ₋₁ F′ₙ₋₁
            iF = inv(Fₙ₋₁)
            Fₙ₋₁ = one(Fₙ₋₁)
            F′ₙ₋₁ *= iF
            lmul!(iF, view(F, i:nF))
            lmul!(iF, view(F′, i:nF))
        end

        F[i] = Fₙ
        F′[i] = F′ₙ
        Fₙ = Fₙ₋₁
        F′ₙ = F′ₙ₋₁
    end
end

function coulomb_upward_recurrence!(G, G′, x⁻¹::T, η, ℓ::UnitRange) where T
    Gₙ = G[1]
    G′ₙ = G′[1]
    η² = η^2
    for i = 2:length(ℓ)
        n = ℓ[i]

        S = (n+1)*x⁻¹ + η/(n+1)
        R = √(1 + η²/(n+1)^2)

        G[i] = Gₙ₊₁ = (S*Gₙ - G′ₙ)/R
        G′[i] = G′ₙ₊₁ = R*Gₙ - S*Gₙ₊₁

        Gₙ = Gₙ₊₁
        G′ₙ = G′ₙ₊₁
    end
end

coulomb_upward_recurrence!(::Nothing, ::Nothing, args...; _...) = nothing

# * Driver

function coulombs!(F::FF, F′::FF, G::GG, G′::GG, x::T, η::T, ℓ::UnitRange; verbosity=0, kwargs...) where {FF,GG,T<:Number}
    verbose = get(kwargs, :verbose, false)
    ρtp = turning_point(η, first(ℓ))
    ρtp > x && verbosity > 0 &&
        @warn "Turning point ρ_TP = $(ρtp) > $(x) which may result in loss of accuracy, consider decreasing first ℓ below $(first(ℓ))."

    if iszero(x)
        F .= false
        for (i,ℓ) in enumerate(ℓ)
            C = coulomb_normalization(η, ℓ)
            F′[i] = (ℓ+1)*C*x^ℓ
            G[i] = 1/((2ℓ+1)*C) # This seems wrong for ℓ ≥ 1
            G′[i] = -T(Inf) # This is not necessarily true
        end
        return
    end

    verbosity > 1 && @show x, η, ℓ
    cf1,n,_,s,converged = coulomb_fraction1(x, η, ℓ[end]; verbosity=verbosity, kwargs...)
    verbosity > 0 && @show cf1
    converged || verbosity > 0 && @warn "Consider increasing ℓmax beyond $(ℓ[end])"

    x⁻¹ = inv(x)
    coulomb_downward_recurrence!(F, F′, x⁻¹, η, ℓ, cf1, s; verbosity=verbosity-2, kwargs...)

    cf2,n,_,s,converged = coulomb_fraction2(x, η, ℓ[1], 1; verbosity=verbosity, kwargs...)
    verbosity > 0 && @show cf2

    p,q = real(cf2),imag(cf2)
    Fₘ = F[1]
    F′ₘ = F′[1]
    Gₘ = (F′ₘ-p*Fₘ)/q
    G′ₘ = p*Gₘ - q*Fₘ

    ω₁ = 1/√(F′ₘ*Gₘ - Fₘ*G′ₘ)
    verbosity > 0 && @show ω₁

    lmul!(ω₁, F)
    lmul!(ω₁, F′)

    if !isnothing(G)
        G[1] = ω₁*Gₘ
        G′[1] = ω₁*G′ₘ

        coulomb_upward_recurrence!(G, G′, x⁻¹, η, ℓ)
    end
end

# * Interface

"""
    coulombs!(F, F′, G, G′, x::AbstractVector, η, ℓ; kwargs...)

Loop through all values of `x` and compute all Coulomb
functions, storing the results in the preallocated matrices `F`, `F′`,
`G`, `G′`.
"""
function coulombs!(F, F′, G, G′, x::AbstractVector, η::Number, ℓ::AbstractVector; kwargs...)
    size(F,1) == size(F′,1) == size(G,1) == size(G′,1) == length(x) &&
        size(F,2) == size(F′,2) == size(G,2) == size(G′,2) == length(ℓ) ||
        throw(DimensionError("The dimension of the output arrays do not agree"))
    for (i,x) in enumerate(x)
        coulombs!(view(F, i, :),
                  view(F′, i, :),
                  view(G, i, :),
                  view(G′, i, :),
                  x, η, ℓ; kwargs...)
    end
end

"""
    coulombs!(F, F′, G, G′, r, Z, k::AbstractVector, ℓ; kwargs...)

Loop through all values of `k` and compute all Coulomb
functions, storing the results in the preallocated matrices `F`, `F′`,
`G`, `G′`. `x=k*r` and `η≡-Z/k`, i.e. `Z>0` for attractive potentials.
"""
function coulombs!(F, F′, G, G′, r::Number, Z::Number, k::AbstractVector, ℓ::AbstractVector; kwargs...)
    size(F,1) == size(F′,1) == size(G,1) == size(G′,1) == length(k) &&
        size(F,2) == size(F′,2) == size(G,2) == size(G′,2) == length(ℓ) ||
        throw(DimensionError("The dimension of the output arrays do not agree"))
    for (i,k) in enumerate(k)
        x = k*r
        η = -Z/k
        coulombs!(view(F, i, :),
                  view(F′, i, :),
                  view(G, i, :),
                  view(G′, i, :),
                  x, η, ℓ; kwargs...)
    end
end

"""
    coulombs(x::AbstractVector, η, ℓ; kwargs...)

Convenience wrapper around [`coulombs!`](@ref) that preallocates output
matrices of the appropriate dimensions.
"""
function coulombs(x::T, η::T, ℓ::AbstractVector; kwargs...) where T
    nℓ = length(ℓ)
    F = zeros(T, nℓ)
    F′ = zeros(T, nℓ)
    G = zeros(T, nℓ)
    G′ = zeros(T, nℓ)
    coulombs!(F, F′, G, G′, x, η, ℓ; kwargs...)
    F, F′, G, G′
end

"""
    coulombs(x::AbstractVector, η, ℓ; kwargs...)

Convenience wrapper around [`coulombs!`](@ref) that preallocates output
matrices of the appropriate dimensions.
"""
function coulombs(x::AbstractVector{T}, η::T, ℓ::AbstractVector; kwargs...) where T
    nx = length(x)
    nℓ = length(ℓ)
    F = zeros(T, nx, nℓ)
    F′ = zeros(T, nx, nℓ)
    G = zeros(T, nx, nℓ)
    G′ = zeros(T, nx, nℓ)
    coulombs!(F, F′, G, G′, x, η, ℓ; kwargs...)
    F, F′, G, G′
end

"""
    coulombs(r, Z, k::AbstractVector, ℓ; kwargs...)

Convenience wrapper around [`coulombs!`](@ref) that preallocates output
matrices of the appropriate dimensions.
"""
function coulombs(r::T, Z::T, k::AbstractVector{T}, ℓ::AbstractVector; kwargs...) where T
    nk = length(k)
    nℓ = length(ℓ)
    F = zeros(T, nk, nℓ)
    F′ = zeros(T, nk, nℓ)
    G = zeros(T, nk, nℓ)
    G′ = zeros(T, nk, nℓ)
    coulombs!(F, F′, G, G′, r, Z, k, ℓ; kwargs...)
    F, F′, G, G′
end

coulombs(x, η, nℓ::Integer; kwargs...) =
    coulombs(x, η, 0:nℓ-1; kwargs...)

coulombs(x, Z, k, nℓ::Integer; kwargs...) =
    coulombs(x, Z, k, 0:nℓ-1; kwargs...)

export coulombs!, coulombs
