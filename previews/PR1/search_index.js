var documenterSearchIndex = {"docs":
[{"location":"spherical_bessel_functions/#Spherical-Bessel-Functions","page":"Spherical Bessel Functions","title":"Spherical Bessel Functions","text":"","category":"section"},{"location":"spherical_bessel_functions/","page":"Spherical Bessel Functions","title":"Spherical Bessel Functions","text":"bessels\nbessels!\nSphericalBesselFunctions.bessel_fraction\nSphericalBesselFunctions.bessel_downward_recurrence!\nSphericalBesselFunctions.neumann_upward_recurrence!\nSphericalBesselFunctions.reflect!\nSphericalBesselFunctions.powneg1","category":"page"},{"location":"spherical_bessel_functions/#SphericalBesselFunctions.bessels","page":"Spherical Bessel Functions","title":"SphericalBesselFunctions.bessels","text":"bessels(x::AbstractVector, nℓ; kwargs...)\n\nConvenience wrapper around bessels! that preallocates output matrices of the appropriate dimensions.\n\n\n\n\n\n","category":"function"},{"location":"spherical_bessel_functions/#SphericalBesselFunctions.bessels!","page":"Spherical Bessel Functions","title":"SphericalBesselFunctions.bessels!","text":"bessels!(j, j′, y, y′, x)\n\nCompute all spherical Bessel (regular: j_n) and Neumann (irregular: y_n) functions and their derivatives at x and store the results in the vectors j, j′, y, y′. If only the Bessel functions or the Neumann functions are of interest, the other pair of arrays can be substituted by nothing. However, it is not possible to compute only the functions but not the derivatives, since they are generated using the following recurrence relations:\n\nbeginaligned\ng_n-1 = fracn+1x g_n + g_n\ng_n-1 = fracn-1x g_n-1 - g_n\nendaligned\n\nThese recurrence relations are employed in a downward fashion for the Bessel functions and an upward fashion for the Neumann functions.\n\nIt is assumed that all passed arrays are of the same lengths (not checked).\n\n\n\n\n\nbessels!(j, j′, y, y′, x::AbstractVector; kwargs...)\n\nLoop through all values of x and compute all Bessel and Neumann functions, storing the results in the preallocated matrices j, j′, y, y′.\n\n\n\n\n\n","category":"function"},{"location":"spherical_bessel_functions/#SphericalBesselFunctions.bessel_fraction","page":"Spherical Bessel Functions","title":"SphericalBesselFunctions.bessel_fraction","text":"bessel_fraction(x, n)\n\nEvaluate the continued fraction for the spherical Bessel function\n\nfracj_n(x)j_n(x) =\nfracnx -\nfrac1T_n+1(x)-frac1T_n+2(x)-frac1T_k(x)-\n\nwhere\n\nT_k(x) = frac2k+1x\n\n\n\n\n\n","category":"function"},{"location":"spherical_bessel_functions/#SphericalBesselFunctions.bessel_downward_recurrence!","page":"Spherical Bessel Functions","title":"SphericalBesselFunctions.bessel_downward_recurrence!","text":"bessel_downward_recurrence!(j, j′, x⁻¹, sinc, cosc, nmax, cf1, s; large)\n\nGiven the logarithmic derivative cf1=j′[end]/j[end] (computed using bessel_fraction), and the sign s, fill in all lower orders using the downward recurrence\n\nbeginaligned\ng_n-1 = S_n+1 g_n + g_n \ng_n-1 = S_n-1 g_n-1 - g_n \nS_n = fracnx\nendaligned\n\nand normalize the whole series using the analytically known expressions\n\ntagDLMF10493\nbeginaligned\nj_0(z) = fracsin zz \nj_0(z) = fraccos zz - fracsin zz^2\nendaligned\n\nIf at any point during the recurrence, j[i]>large, renormalize j′[i:end] ./= j[i] and j[i:end] ./= j[i] to avoid overflow, and the continue the recurrence. This is very useful for small x and large n.\n\nnmax allows starting the recurrence for a higher n than for which there is space allocated in j and j′ (this can help convergence for the terms which are actually of interest).\n\n\n\n\n\n","category":"function"},{"location":"spherical_bessel_functions/#SphericalBesselFunctions.neumann_upward_recurrence!","page":"Spherical Bessel Functions","title":"SphericalBesselFunctions.neumann_upward_recurrence!","text":"neumann_upward_recurrence(y, y′, x⁻¹, sinc, cosc)\n\nGenerate the irregular Neumann functions using (stable) upward recurrence\n\nbeginaligned\ng_n+1 = S_n g_n - g_n \ng_n+1 = g_n - S_n+2 g_n+1\nendaligned\n\nstarting from the analytically known expressions\n\ntagDLMF10495\nbeginaligned\ny_0(z) = -fraccos zz \ny_0(z) = fracsin zz + fraccos zz^2\nendaligned\n\n\n\n\n\n","category":"function"},{"location":"spherical_bessel_functions/#SphericalBesselFunctions.reflect!","page":"Spherical Bessel Functions","title":"SphericalBesselFunctions.reflect!","text":"reflect!(g, g′, offset=0)\n\nAffect the reflection symmetries\n\ntagDLMF104714\nbeginaligned\nj_n(-z) = (-)^n j_n(z) \ny_n(-z) = (-)^n+1y_n(z) \nj_n(-z) = (-)^n-1 j_n(z) \ny_n(-z) = (-)^n y_n(z)\nendaligned\n\nwhere the +1 can be added using offset.\n\n\n\n\n\n","category":"function"},{"location":"spherical_bessel_functions/#SphericalBesselFunctions.powneg1","page":"Spherical Bessel Functions","title":"SphericalBesselFunctions.powneg1","text":"powneg1(m)\n\nReturns an integer power of negative unity, i.e. (-)^n.\n\n\n\n\n\n","category":"function"},{"location":"continued_fractions/#Continued-Fractions","page":"Continued Fractions","title":"Continued Fractions","text":"","category":"section"},{"location":"continued_fractions/","page":"Continued Fractions","title":"Continued Fractions","text":"Sometimes it is difficult to know a priori how many iterations are needed for the continued fractions to converge. In the figure below, we illustrate how many iterations were needed to reach convergence for SphericalBesselFunctions.coulomb_fraction1 (solid lines) and SphericalBesselFunctions.coulomb_fraction2 (dashed lines), as a function of x, for a range of values of eta (pm10^n, n=-13) and lambda=0. Shown in black is the number of iterations required to converge SphericalBesselFunctions.bessel_fraction. This plot should be compared with figure 1 of","category":"page"},{"location":"continued_fractions/","page":"Continued Fractions","title":"Continued Fractions","text":"Barnett, A. (1982). Continued-fraction evaluation of Coulomb functions F_lambda(eta x), G_lambda(eta x) and their derivatives. Journal of Computational Physics, 46(2), 171–188. 10.1016/0021-9991(82)90012-2","category":"page"},{"location":"continued_fractions/","page":"Continued Fractions","title":"Continued Fractions","text":"(Image: Continued fractions)","category":"page"},{"location":"continued_fractions/","page":"Continued Fractions","title":"Continued Fractions","text":"SphericalBesselFunctions.lentz_thompson\nSphericalBesselFunctions.steed_kahan","category":"page"},{"location":"continued_fractions/#SphericalBesselFunctions.lentz_thompson","page":"Continued Fractions","title":"SphericalBesselFunctions.lentz_thompson","text":"lentz_thompson(b₀, a, b; max_iter=20_000)\n\nLentz–Thompson algorithm for the forward evaluation of continued fractions:\n\nf_n = b_0 + fraca_1b_1+fraca_2b_2+fraca_n-1b_n-1+fraca_nb_n\n\nAs described in\n\nBarnett, A. R. (1996). The Calculation of Spherical Bessel and Coulomb Functions. In (Eds.), Computational Atomic Physics (pp. 181–202). Springer Berlin Heidelberg.\n\nThe number of iterations required is roughly proportional to the argument of the Bessel, according to Barnett (1996). For the Coulomb functions, the number of iterations required is a bit more involved, see\n\nBarnett, A. (1982). Continued-fraction evaluation of Coulomb functions F_lambda(eta x), G_lambda(eta x) and their derivatives. Journal of Computational Physics, 46(2), 171–188. 10.1016/0021-9991(82)90012-2\n\n\n\n\n\n","category":"function"},{"location":"continued_fractions/#SphericalBesselFunctions.steed_kahan","page":"Continued Fractions","title":"SphericalBesselFunctions.steed_kahan","text":"steed_kahan(b₀, a, b)\n\nSteed's algorithm for the forward evaluation of continued fractions, with Kahan summation to avoid loss of accuracy.\n\nf_n = b_0 + fraca_1b_1+fraca_2b_2+fraca_n-1b_n-1+fraca_nb_n\n\nAs described in\n\nThompson, I., & Barnett, A. (1986). Coulomb and Bessel functions of complex arguments and order. Journal of Computational Physics, 64(2), 490–509. 10.1016/0021-9991(86)90046-x\n\n\n\n\n\n","category":"function"},{"location":"coulomb_functions/#Coulomb-Functions","page":"Coulomb Functions","title":"Coulomb Functions","text":"","category":"section"},{"location":"coulomb_functions/","page":"Coulomb Functions","title":"Coulomb Functions","text":"coulombs\ncoulombFs\ncoulombs!\nSphericalBesselFunctions.coulomb_fraction1\nSphericalBesselFunctions.coulomb_fraction2\nSphericalBesselFunctions.coulomb_downward_recurrence!\nSphericalBesselFunctions.coulomb_upward_recurrence!\nSphericalBesselFunctions.turning_point\nSphericalBesselFunctions.coulomb_normalization","category":"page"},{"location":"coulomb_functions/#SphericalBesselFunctions.coulombs","page":"Coulomb Functions","title":"SphericalBesselFunctions.coulombs","text":"coulombs(x::AbstractVector, η, ℓ; kwargs...)\n\nConvenience wrapper around coulombs! that preallocates output matrices of the appropriate dimensions.\n\n\n\n\n\ncoulombs(x::AbstractVector, η, ℓ; kwargs...)\n\nConvenience wrapper around coulombs! that preallocates output matrices of the appropriate dimensions.\n\n\n\n\n\ncoulombs(r, Z, k::AbstractVector, ℓ; kwargs...)\n\nConvenience wrapper around coulombs! that preallocates output matrices of the appropriate dimensions.\n\n\n\n\n\n","category":"function"},{"location":"coulomb_functions/#SphericalBesselFunctions.coulombFs","page":"Coulomb Functions","title":"SphericalBesselFunctions.coulombFs","text":"coulombFs(x::AbstractVector, η, ℓ; kwargs...)\n\nConvenience wrapper around coulombs! that preallocates output matrices of the appropriate dimensions, only computing the regular functions F_ell(etax) and F_ell(etax).\n\n\n\n\n\n","category":"function"},{"location":"coulomb_functions/#SphericalBesselFunctions.coulombs!","page":"Coulomb Functions","title":"SphericalBesselFunctions.coulombs!","text":"coulombs!(F, F′, G, G′, x, η, ℓ)\n\nMain routine that computes the continued fractions and the recurrence relations to generate the regular functions F_ll(tax) and F_ll(tax) and the irregular functions G_ll(tax) and G_ll(tax) (computation of the latter can be elided by passing nothing for G and G′).\n\n\n\n\n\ncoulombs!(F, F′, G, G′, x::AbstractVector, η, ℓ; kwargs...)\n\nLoop through all values of x and compute all Coulomb functions, storing the results in the preallocated matrices F, F′, G, G′.\n\n\n\n\n\ncoulombs!(F, F′, G, G′, r, Z, k::AbstractVector, ℓ; kwargs...)\n\nLoop through all values of k and compute all Coulomb functions, storing the results in the preallocated matrices F, F′, G, G′. x=k*r and η≡-Z/k, i.e. Z>0 for attractive potentials.\n\n\n\n\n\n","category":"function"},{"location":"coulomb_functions/#SphericalBesselFunctions.coulomb_fraction1","page":"Coulomb Functions","title":"SphericalBesselFunctions.coulomb_fraction1","text":"coulomb_fraction1(x, η, n)\n\nEvaluate the continued fraction for the logarithmic derivative of the Coulomb function\n\nfracF_n(etax)F_n(etax) =\nS_n+1(etax) -\nfracR_n+1^2(eta)T_n+1(etax)-fracR_n+2^2(eta)T_n+2(etax)-fracR_k^2(eta)T_k(etax)-\n\nwhere\n\nbeginaligned\nR_k(eta) = sqrt1+fraceta^2k^2 qquad\nS_k(etax) = frackx + fracetak \nT_k(etax) = S_k(etax) + S_k+1(etax) = (2k+1)left(frac1x + fracetak^2+kright)\nendaligned\n\n\n\n\n\n","category":"function"},{"location":"coulomb_functions/#SphericalBesselFunctions.coulomb_fraction2","page":"Coulomb Functions","title":"SphericalBesselFunctions.coulomb_fraction2","text":"coulomb_fraction2(x, η, n, ω)\n\nEvaluate the continued fraction for the logarithmic derivative of the Hankel function H^omega_lambda(etax) = G_lambda(etax) + mathrmiomega F_lambda(etax)\n\nfracH^omega_n(etax)H^omega_n(etax) =\np+mathrmiq =\nmathrmiomegaleft(1-fracetaxright) +\nfracmathrmiomegax\nfracac2(x-eta+mathrmiomega)+\nfrac(a+1)(c+1)2(x-eta+2mathrmiomega)+\n\nwhere\n\nbeginaligned\na = 1+n+mathrmiomegaeta \nb = 2n + 2 \nc = -n + mathrmiomegaeta\nendaligned\n\n\n\n\n\n","category":"function"},{"location":"coulomb_functions/#SphericalBesselFunctions.coulomb_downward_recurrence!","page":"Coulomb Functions","title":"SphericalBesselFunctions.coulomb_downward_recurrence!","text":"coulomb_downward_recurrence!(F, F′, x⁻¹, η, ℓ, cf1, s; large)\n\nGiven the logarithmic derivative cf1=F′[end]/F[end] (computed using coulomb_fraction1), and the sign s, fill in all lower orders using the downward recurrence\n\nbeginaligned\nw_n-1 =\n(S_n w_n + w_n)R_n \nw_n-1 =\nS_nw_n-1 - R_n w_n \nR_k = sqrt1 + fraceta^2k^2 \nS_k = frackx + fracetak\nendaligned\n\nIf at any point during the recurrence, F[i]>large, renormalize F′[i:end] ./= F[i] and F[i:end] ./= F[i] to avoid overflow, and the continue the recurrence. This is very useful for small x and large ell.\n\n\n\n\n\n","category":"function"},{"location":"coulomb_functions/#SphericalBesselFunctions.coulomb_upward_recurrence!","page":"Coulomb Functions","title":"SphericalBesselFunctions.coulomb_upward_recurrence!","text":"coulomb_upward_recurrence(G, G′, x⁻¹, η, ℓ)\n\nGenerate the irregular Coulomb functions using (stable) upward recurrence\n\nbeginaligned\nw_n+1 = (S_n+1 w_n - w_n)R_n+1 \nw_n+1 = R_n+1w_n - S_n+1 w_n+1\nendaligned\n\nstarting from G[1] and G′[1], which we find using the Wronskian\n\nF_lambda(etaz)G_lambda(etaz) - F_lambda(etaz)G_lambda(etaz) = 1\n\nand the logarithmic derivative of G_lambda(etaz) + mathrmiF_lambda(etaz) (see coulomb_fraction2).\n\n\n\n\n\n","category":"function"},{"location":"coulomb_functions/#SphericalBesselFunctions.turning_point","page":"Coulomb Functions","title":"SphericalBesselFunctions.turning_point","text":"turning_point(η, ℓ)\n\nComputes the turning point rho_mathrmTP of the Coulomb functions, which separates the monotonic region rhorho_mathrmTP from the oscillatory region rhorho_mathrmTP.\n\nSee Eq. (2) of\n\nBarnett, A. (1982). COULFG: Coulomb and Bessel Functions and Their Derivatives, for Real Arguments, by Steed's Method. Computer Physics Communications, 27(2), 147–166. 10.1016/0010-4655(82)90070-4\n\nor\n\nBarnett, A. R. (1996). The Calculation of Spherical Bessel and Coulomb Functions. In (Eds.), Computational Atomic Physics (pp. 181–202). : Springer Berlin Heidelberg.\n\n\n\n\n\n","category":"function"},{"location":"coulomb_functions/#SphericalBesselFunctions.coulomb_normalization","page":"Coulomb Functions","title":"SphericalBesselFunctions.coulomb_normalization","text":"coulomb_normalization(η, ℓ)\n\nCompute the Coulombic normalization constant\n\ntagDLMF3325\nC_ell(eta) =\nfrac2^ell mathrme^-pieta2Gamma(ell+1+mathrmieta)Gamma(2ell+2)\n\n\n\n\n\n","category":"function"},{"location":"#SphericalBesselFunctions.jl","page":"Home","title":"SphericalBesselFunctions.jl","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"This library provide efficient and accurate implementations of the regular and irregular spherical Bessel functions [j_n(z) and y_n(z)] and Coulomb functions [F_lambda(etaz) and G_lambda(etaz)], and their derivatives. The former are solutions to the radial part of the Helmholtz equation in spherical coordinates:","category":"page"},{"location":"","page":"Home","title":"Home","text":"tagDLMF10471\nz^2w + 2zw + z^2 - n(n+1) = 0","category":"page"},{"location":"","page":"Home","title":"Home","text":"whereas the latter obey","category":"page"},{"location":"","page":"Home","title":"Home","text":"tagDLMF3321\nw + left\n1 - frac2etaz -\nfraclambda(lambda+1)z^2\nrightw = 0","category":"page"},{"location":"","page":"Home","title":"Home","text":"The spherical Bessel functions are related to the Coulomb functions as","category":"page"},{"location":"","page":"Home","title":"Home","text":"tagDLMF3353\nbeginaligned\nj_n(z) = fracF_n(0z)z \ny_n(z) = -fracG_n(0z)z\nendaligned","category":"page"},{"location":"#Usage","page":"Home","title":"Usage","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"x = range(-12, stop=12, length=1000)\nnℓ = 10\nj, j′, y, y′ = bessels(x, nℓ)","category":"page"},{"location":"","page":"Home","title":"Home","text":"(Image: Simple Bessel example)","category":"page"},{"location":"","page":"Home","title":"Home","text":"x = range(0, stop=12, length=1000)\nnℓ = 10\nη = -1.0\nF, F′, G, G′ = coulombs(x, η, nℓ)","category":"page"},{"location":"","page":"Home","title":"Home","text":"(Image: Simple Coulomb example)","category":"page"},{"location":"#Accuracy","page":"Home","title":"Accuracy","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"To check the accuracy, we compare with SpecialFunctions.jl, which however does not provide the spherical Bessel functions but the ordinary (cylindrical) ones. They are however related as","category":"page"},{"location":"","page":"Home","title":"Home","text":"tagDLMF104734\nbeginaligned\nj_n(z) equiv sqrtfracpi2z\nJ_n+frac12(z)\ny_n(z) equiv sqrtfracpi2z\nY_n+frac12(z)\nendaligned","category":"page"},{"location":"","page":"Home","title":"Home","text":"Furthermore, SpecialFunctions.jl does not provide the derivatives out-of-the-box, but they are easily found using the recurrence relation.:","category":"page"},{"location":"","page":"Home","title":"Home","text":"f_n(z) = f_n-1(z) - fracn+1z f_n(z)\nqquad\nf_n = j_ny_n","category":"page"},{"location":"","page":"Home","title":"Home","text":"This time, we investigate a larger domain of parameters, but avoid smaller values of x than 01 since that does not seem to work in SpecialFunctions.jl (SphericalBesselFunctions.jl works at xleq0 as well):","category":"page"},{"location":"","page":"Home","title":"Home","text":"nx = 1001\nx = 10 .^ range(-1, stop=4, length=nx)\nnℓ = 105\nj, j′, y, y′ = bessels(x, nℓ)","category":"page"},{"location":"","page":"Home","title":"Home","text":"(Image: Accuracy)","category":"page"},{"location":"","page":"Home","title":"Home","text":"We note that SphericalBesselFunctions.jl seems to be ~3 times faster than SpecialFunctions.jl when evaluating all Bessel functions for a fixed value of x, most likely due to extra processing taking place when SpecialFunctions.jl does not compute all values simultaneously, but one order at a time:","category":"page"},{"location":"","page":"Home","title":"Home","text":"SphericalBesselFunctions.jl:","category":"page"},{"location":"","page":"Home","title":"Home","text":"BenchmarkTools.Trial:\n  memory estimate:  0 bytes\n  allocs estimate:  0\n  --------------\n  minimum time:     76.331 μs (0.00% GC)\n  median time:      77.057 μs (0.00% GC)\n  mean time:        80.337 μs (0.00% GC)\n  maximum time:     197.715 μs (0.00% GC)\n  --------------\n  samples:          10000\n  evals/sample:     1","category":"page"},{"location":"","page":"Home","title":"Home","text":"SpecialFunctions.jl:","category":"page"},{"location":"","page":"Home","title":"Home","text":"BenchmarkTools.Trial:\n  memory estimate:  0 bytes\n  allocs estimate:  0\n  --------------\n  minimum time:     170.295 μs (0.00% GC)\n  median time:      185.132 μs (0.00% GC)\n  mean time:        193.312 μs (0.00% GC)\n  maximum time:     442.129 μs (0.00% GC)\n  --------------\n  samples:          10000\n  evals/sample:     1","category":"page"},{"location":"","page":"Home","title":"Home","text":"Finally, the agreement for the (irregular) Neumann functions for x  100 is terrible, unclear why. But they are irregular (diverging) as x tends to zero.","category":"page"}]
}
