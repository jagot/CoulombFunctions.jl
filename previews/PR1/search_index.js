var documenterSearchIndex = {"docs":
[{"location":"#SphericalBesselFunctions.jl","page":"Home","title":"SphericalBesselFunctions.jl","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"#Usage","page":"Home","title":"Usage","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"x = range(-12, stop=12, length=1000)\nnℓ = 10\nj, j′, y, y′ = bessels(x, nℓ)","category":"page"},{"location":"","page":"Home","title":"Home","text":"(Image: Simple example)","category":"page"},{"location":"#Accuracy","page":"Home","title":"Accuracy","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"To check the accuracy, we compare with SpecialFunctions.jl, which however does not provide the spherical Bessel functions but the ordinary (cylindrical) ones. They are however related as","category":"page"},{"location":"","page":"Home","title":"Home","text":"beginaligned\nj_n(z) equiv sqrtfracpi2z\nJ_n+frac12(z)\ny_n(z) equiv sqrtfracpi2z\nY_n+frac12(z)\nendaligned","category":"page"},{"location":"","page":"Home","title":"Home","text":"Furthermore, SpecialFunctions.jl does not provide the derivatives out-of-the-box, but they are easily found using the recurrence relation.:","category":"page"},{"location":"","page":"Home","title":"Home","text":"f_n(z) = f_n-1(z) - fracn+1z f_n(z)\nqquad\nf_n = j_ny_n","category":"page"},{"location":"","page":"Home","title":"Home","text":"This time, we investigate a larger domain of parameters, but avoid smaller values of x than 01 since that does not seem to work in SpecialFunctions.jl (SphericalBesselFunctions.jl works at xleq0 as well; however the combination large n and small x is still problematic):","category":"page"},{"location":"","page":"Home","title":"Home","text":"nx = 1001\nx = 10 .^ range(-1, stop=4, length=nx)\nnℓ = 105\nj, j′, y, y′ = bessels(x, nℓ)","category":"page"},{"location":"","page":"Home","title":"Home","text":"(Image: Accuracy)","category":"page"},{"location":"","page":"Home","title":"Home","text":"We note that SphericalBesselFunctions.jl seems to be ~3 times faster than SpecialFunctions.jl when evaluating all Bessel functions for a fixed value of x, most likely due to extra processing taking place when SpecialFunctions.jl does not compute all values simultaneously, but one order at a time:","category":"page"},{"location":"","page":"Home","title":"Home","text":"SphericalBesselFunctions.jl:","category":"page"},{"location":"","page":"Home","title":"Home","text":"BenchmarkTools.Trial:\n  memory estimate:  0 bytes\n  allocs estimate:  0\n  --------------\n  minimum time:     76.331 μs (0.00% GC)\n  median time:      77.057 μs (0.00% GC)\n  mean time:        80.337 μs (0.00% GC)\n  maximum time:     197.715 μs (0.00% GC)\n  --------------\n  samples:          10000\n  evals/sample:     1","category":"page"},{"location":"","page":"Home","title":"Home","text":"SpecialFunctions.jl:","category":"page"},{"location":"","page":"Home","title":"Home","text":"BenchmarkTools.Trial:\n  memory estimate:  0 bytes\n  allocs estimate:  0\n  --------------\n  minimum time:     170.295 μs (0.00% GC)\n  median time:      185.132 μs (0.00% GC)\n  mean time:        193.312 μs (0.00% GC)\n  maximum time:     442.129 μs (0.00% GC)\n  --------------\n  samples:          10000\n  evals/sample:     1","category":"page"},{"location":"","page":"Home","title":"Home","text":"Finally, the agreement for the (irregular) Neumann functions for x  100 is terrible, unclear why. But they are irregular (diverging) as x tends to zero.","category":"page"},{"location":"#Reference","page":"Home","title":"Reference","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Modules = [SphericalBesselFunctions]","category":"page"},{"location":"#SphericalBesselFunctions.bessel_fraction-Union{Tuple{T}, Tuple{T, Integer}} where T","page":"Home","title":"SphericalBesselFunctions.bessel_fraction","text":"bessel_fraction(x, n)\n\nEvaluated the continued fraction for the spherical Bessel function\n\nfracj_n(x)j_n(x) =\nfracnx -\nfrac1T_n+1(x)-frac1T_n+2(x)-frac1T_k-\n\nwhere\n\nT_k(x) = frac2k+1x\n\n\n\n\n\n","category":"method"},{"location":"#SphericalBesselFunctions.bessels!-Tuple{Any, Any, Any, Any, AbstractVector}","page":"Home","title":"SphericalBesselFunctions.bessels!","text":"bessels!(j, j′, y, y′, x::AbstractVector; kwargs...)\n\nLoop through all values of x and compute all Bessel and Neumann functions, storing the results in the preallocated matrices j, j′, y, y′.\n\n\n\n\n\n","category":"method"},{"location":"#SphericalBesselFunctions.bessels!-Union{Tuple{T}, Tuple{Y}, Tuple{J}, Tuple{J, J, Y, Y, T}} where {J, Y, T<:Number}","page":"Home","title":"SphericalBesselFunctions.bessels!","text":"bessels!(j, j′, y, y′, x)\n\nCompute all spherical Bessel (regular: j_n) and Neumann (irregular: y_n) functions and their derivatives at x and store the results in the vectors j, j′, y, y′. If only the Bessel functions or the Neumann functions are of interest, the other pair of arrays can be substituted by nothing. However, it is not possible to compute only the functions but not the derivatives, since they are generated using the following recurrence relations:\n\nbeginaligned\ng_n-1 = fracn+1x g_n + g_n\ng_n-1 = fracn-1x g_n-1 - g_n\nendaligned\n\nThese recurrence relations are employed in a downward fashion for the Bessel functions and an upward fashion for the Neumann functions.\n\nIt is assumed that all passed arrays are of the same lengths (not checked).\n\n\n\n\n\n","category":"method"},{"location":"#SphericalBesselFunctions.bessels-Union{Tuple{T}, Tuple{AbstractVector{T}, Any}} where T","page":"Home","title":"SphericalBesselFunctions.bessels","text":"bessels(x::AbstractVector, nℓ; kwargs...)\n\nConvenience wrapper around bessels! that preallocates output matrices of the appropriate dimensions.\n\n\n\n\n\n","category":"method"},{"location":"#SphericalBesselFunctions.coulombs!-Tuple{Any, Any, Any, Any, AbstractVector, Any, AbstractVector}","page":"Home","title":"SphericalBesselFunctions.coulombs!","text":"coulombs!(F, F′, G, G′, x::AbstractVector; kwargs...)\n\nLoop through all values of x and compute all Coulomb functions, storing the results in the preallocated matrices F, F′, G, G′.\n\n\n\n\n\n","category":"method"},{"location":"#SphericalBesselFunctions.coulombs-Union{Tuple{T}, Tuple{AbstractVector{T}, Any, AbstractVector}} where T","page":"Home","title":"SphericalBesselFunctions.coulombs","text":"coulombs(x::AbstractVector, nℓ; kwargs...)\n\nConvenience wrapper around coulombs! that preallocates output matrices of the appropriate dimensions.\n\n\n\n\n\n","category":"method"},{"location":"#SphericalBesselFunctions.lentz_thompson-Union{Tuple{T}, Tuple{T, Any, Any}} where T","page":"Home","title":"SphericalBesselFunctions.lentz_thompson","text":"lentz_thompson(b₀, a, b; max_iter=20_000)\n\nLentz–Thompson algorithm for the forward evaluation of continued fractions:\n\nf_n = b_0 + fraca_1b_1+fraca_2b_2+fraca_n-1b_n-1+fraca_nb_n\n\nAs described in\n\nBarnett, A. R. (1996). The Calculation of Spherical Bessel and Coulomb Functions. In (Eds.), Computational Atomic Physics (pp. 181–202). Springer Berlin Heidelberg.\n\nThe number of iterations required is roughly proportional to the argument of the Bessel/Coulomb function, according to Barnett (1996) and\n\nBarnett, A. (1982). Continued-fraction evaluation of coulomb functions F_lambda(eta x), Glambda(eta x) and their derivatives. Journal of Computational Physics, 46(2), 171–188. http://dx.doi.org/10.1016/0021-9991(82)90012-2\n\n\n\n\n\n","category":"method"},{"location":"#SphericalBesselFunctions.powneg1-Tuple{Integer}","page":"Home","title":"SphericalBesselFunctions.powneg1","text":"powneg1(m)\n\nReturns an integer power of negative unity, i.e. (-)^n.\n\n\n\n\n\n","category":"method"},{"location":"#SphericalBesselFunctions.reflect!","page":"Home","title":"SphericalBesselFunctions.reflect!","text":"reflect!(g, g′, offset=0)\n\nAffect the reflection symmetries\n\nj_n(-z) = (-)^n j_n(z) qquad\ny_n(-z) = (-)^n+1y_n(z)\n\nwhere the +1 can be added using offset.\n\n\n\n\n\n","category":"function"},{"location":"#SphericalBesselFunctions.steed_kahan-Union{Tuple{T}, Tuple{T, Any, Any}} where T","page":"Home","title":"SphericalBesselFunctions.steed_kahan","text":"steed_kahan(b₀, a, b)\n\nSteed's algorithm for the forward evaluation of continued fractions, with Kahan summation to avoid loss of accuracy.\n\nf_n = b_0 + fraca_1b_1+fraca_2b_2+fraca_n-1b_n-1+fraca_nb_n\n\nAs described in\n\nThompson, I., & Barnett, A. (1986). Coulomb and Bessel functions of complex arguments and order. Journal of Computational Physics, 64(2), 490–509. http://dx.doi.org/10.1016/0021-9991(86)90046-x\n\n\n\n\n\n","category":"method"},{"location":"#SphericalBesselFunctions.turning_point-Tuple{Any, Any}","page":"Home","title":"SphericalBesselFunctions.turning_point","text":"turning_point(η, ℓ)\n\nComputes the turning point rho_mathrmTP of the Coulomb functions, which separates the monotonic region rhorho_mathrmTP from the oscillatory region rhorho_mathrmTP.\n\nSee Eq. (2) of\n\nBarnett, A. (1982). COULFG: Coulomb and Bessel Functions and Their Derivatives, for Real Arguments, by Steed's Method. Computer Physics Communications, 27(2), 147–166. http://dx.doi.org/10.1016/0010-4655(82)90070-4\n\nor\n\nBarnett, A. R. (1996). The Calculation of Spherical Bessel and Coulomb Functions. In (Eds.), Computational Atomic Physics (pp. 181–202). : Springer Berlin Heidelberg.\n\n\n\n\n\n","category":"method"}]
}
