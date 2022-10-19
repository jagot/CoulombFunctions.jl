# Continued Fractions

Sometimes it is difficult to know _a priori_ how many iterations are
needed for the continued fractions to converge. In the figure below,
we illustrate how many iterations were needed to reach convergence for
[`SphericalBesselFunctions.coulomb_fraction1`](@ref) (solid lines) and
[`SphericalBesselFunctions.coulomb_fraction2`](@ref) (dashed lines),
as a function of ``x``, for a range of values of ``\eta``
(``\pm10^n``, ``n=-1..3``) and ``\lambda=0``. Shown in black is the
number of iterations required to converge
[`SphericalBesselFunctions.bessel_fraction`](@ref). This plot should
be compared with figure 1 of

- Barnett, A. (1982). Continued-fraction evaluation of Coulomb
  functions ``F_\lambda(\eta, x)``, ``G_\lambda(\eta, x)`` and their
  derivatives. Journal of Computational Physics, 46(2),
  171–188. [10.1016/0021-9991(82)90012-2](http://dx.doi.org/10.1016/0021-9991(82)90012-2)

![Continued fractions](figures/continued-fractions.svg)

```@docs
SphericalBesselFunctions.lentz_thompson
SphericalBesselFunctions.steed_kahan
```
