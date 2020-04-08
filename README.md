# SphericalBesselFunctions

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://jagot.github.io/SphericalBesselFunctions.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://jagot.github.io/SphericalBesselFunctions.jl/dev)
[![Build Status](https://travis-ci.com/jagot/SphericalBesselFunctions.jl.svg?branch=master)](https://travis-ci.com/jagot/SphericalBesselFunctions.jl)
[![Build Status](https://ci.appveyor.com/api/projects/status/github/jagot/SphericalBesselFunctions.jl?svg=true)](https://ci.appveyor.com/project/jagot/SphericalBesselFunctions-jl)
[![Codecov](https://codecov.io/gh/jagot/SphericalBesselFunctions.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/jagot/SphericalBesselFunctions.jl)

A small Julia package for the evaluation of [spherical Bessel
functions](https://dlmf.nist.gov/10#PT4) using the approach described
in

- Barnett, A. R. (1996). The Calculation of Spherical Bessel and
  Coulomb Functions. In (Eds.), Computational Atomic Physics
  (pp. 181–202). Springer Berlin Heidelberg.

The article and accompanying codes (public domain) can be found at
http://www.fresco.org.uk/programs/barnett/index.htm

A fix in the rescaling step occurring in the vicinity of roots of the
`sinc` function was implemented similarly to the suggestion here:
https://github.com/manoharan-lab/holopy/pull/273/files
