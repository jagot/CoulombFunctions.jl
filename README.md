# SphericalBesselFunctions

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://jagot.github.io/SphericalBesselFunctions.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://jagot.github.io/SphericalBesselFunctions.jl/dev/)
[![Build Status](https://github.com/jagot/SphericalBesselFunctions.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/jagot/SphericalBesselFunctions.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/jagot/SphericalBesselFunctions.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/jagot/SphericalBesselFunctions.jl)

A small Julia package for the evaluation of [spherical Bessel
functions](https://dlmf.nist.gov/10#PT4) and [Coulomb
functions](https://dlmf.nist.gov/33#PT2) using the approach described
in

- Barnett, A. R. (1996). The Calculation of Spherical Bessel and
  Coulomb Functions. In (Eds.), Computational Atomic Physics
  (pp. 181–202). Springer Berlin Heidelberg.

The article and accompanying codes (public domain) can be found at
http://www.fresco.org.uk/programs/barnett/index.htm

A fix in the rescaling step occurring in the vicinity of roots of the
`sinc` function was implemented similarly to the suggestion here:
https://github.com/manoharan-lab/holopy/pull/273/files

Alternatives:
  * [Bessels.jl](https://github.com/JuliaMath/Bessels.jl/)
  * [SpecialFunctions.jl](https://github.com/JuliaMath/SpecialFunctions.jl)
