using PyPlot
using Jagot.plotting
plot_style("ggplot")

using PyCall
Cycler = pyimport("cycler")
plt.rc("axes", prop_cycle=plt.rcParams["axes.prop_cycle"])

using CoulombFunctions
import CoulombFunctions: bessel_fraction, coulomb_fraction1, coulomb_fraction2

include(joinpath(@__DIR__, "../test/reference.jl"))
FIGDIR = joinpath(@__DIR__, "src", "figures")

function simple_example()
    nx = 1000
    x = range(-12, stop=12, length=nx)

    nℓ = 10
    j, j′, y, y′ = bessels(x, nℓ)

    cfigure("simple example",figsize=(7,6)) do
        csubplot(221,nox=true) do
            plot(x, j)
            margins(0, 0.1)
            ylim(-0.5, 0.5)
            ylabel(L"j_n(x)")
        end
        csubplot(222,nox=true) do
            plot(x, y)
            margins(0, 0.1)
            ylim(-0.5, 0.5)
            ylabel(L"y_n(x)")
            axes_labels_opposite(:y)
        end
        csubplot(223) do
            plot(x, j′)
            margins(0, 0.1)
            ylim(-0.5, 0.5)
            xlabel(L"x")
            ylabel(L"j'_n(x)")
        end
        csubplot(224) do
            plot(x, y′)
            margins(0, 0.1)
            ylim(-0.5, 0.5)
            xlabel(L"x")
            ylabel(L"y'_n(x)")
            axes_labels_opposite(:y)
        end
    end
    savefig(joinpath(FIGDIR, "simple-example.svg"))
end

function simple_coulomb_example()
    nx = 1000
    x = range(0.0, stop=20, length=nx)

    nℓ = 10
    η = -1.0
    F, F′, G, G′ = coulombs(x, η, nℓ)

    cfigure("simple coulomb example",figsize=(7,6)) do
        csubplot(221,nox=true) do
            plot(x, F)
            margins(0.05, 0.1)
            ylim(-2.0, 2.0)
            ylabel(L"F_\lambda(\eta,x)")
        end
        csubplot(222,nox=true) do
            plot(x, G)
            margins(0.05, 0.1)
            ylim(-2.0, 2.0)
            ylabel(L"G_\lambda(\eta,x)")
            axes_labels_opposite(:y)
        end
        csubplot(223) do
            plot(x, F′)
            margins(0.05, 0.1)
            ylim(-2.0, 2.0)
            xlabel(L"x")
            ylabel(L"F'_\lambda(\eta,x)")
        end
        csubplot(224) do
            plot(x, G′)
            margins(0.05, 0.1)
            ylim(-2.0, 2.0)
            xlabel(L"x")
            ylabel(L"G'_\lambda(\eta,x)")
            axes_labels_opposite(:y)
        end
    end
    savefig(joinpath(FIGDIR, "simple-coulomb-example.svg"))
end

function accuracy()
    nx = 1001
    x = 10 .^ range(-1, stop=4, length=nx)
    nℓ = 105
    j, j′, y, y′ = @time bessels(x, nℓ)

    sel = 1:min(nℓ,104)

    nℓref = min(nℓ,length(sel))
    jref = zeros(nx, nℓref)
    j′ref = zeros(nx, nℓref)
    yref = zeros(nx, nℓref)
    y′ref = zeros(nx, nℓref)

    @time for (i,x) in enumerate(x)
        bessels_ref!(view(jref, i, :),
                     view(j′ref, i, :),
                     view(yref, i, :),
                     view(y′ref, i, :),
                     x)
    end

    j = j[:,sel]
    j′ = j′[:,sel]
    y = y[:,sel]
    y′ = y′[:,sel]
    jref = jref[:,sel]
    j′ref = j′ref[:,sel]
    yref = yref[:,sel]
    y′ref = y′ref[:,sel]

    cfigure("accuracy", figsize=(9,10)) do
        m,n = 4,3
        ks = reshape(1:m*n, n, m)'

        csubplot(m,n,ks[1,1],nox=true) do
            semilogx(x, j, rasterized=true)
            ylim(-0.5,0.5)
            ylabel(L"j_n(x)")
            title("CoulombFunctions.jl")
        end
        csubplot(m,n,ks[2,1],nox=true) do
            semilogx(x, j′, rasterized=true)
            ylim(-0.5,0.5)
            ylabel(L"j_n'(x)")
        end
        csubplot(m,n,ks[3,1],nox=true) do
            semilogx(x, y, rasterized=true)
            ylim(-0.5,0.5)
            ylabel(L"y_n(x)")
        end
        csubplot(m,n,ks[4,1]) do
            semilogx(x, y′, rasterized=true)
            ylim(-0.5,0.5)
            xlabel(L"x")
            ylabel(L"y_n'(x)")
        end

        csubplot(m,n,ks[1,2],nox=true,noy=true) do
            semilogx(x, jref, rasterized=true)
            ylim(-0.5,0.5)
            title("SpecialFunctions.jl")
        end
        csubplot(m,n,ks[2,2],nox=true,noy=true) do
            semilogx(x, j′ref, rasterized=true)
            ylim(-0.5,0.5)
        end
        csubplot(m,n,ks[3,2],nox=true,noy=true) do
            semilogx(x, yref, rasterized=true)
            ylim(-0.5,0.5)
        end
        csubplot(m,n,ks[4,2],noy=true) do
            semilogx(x, y′ref, rasterized=true)
            ylim(-0.5,0.5)
            xlabel(L"x")
        end

        csubplot(m,n,ks[1,3],nox=true) do
            jerror = abs.(j-jref)
            loglog(x, jerror, rasterized=true)
            ylim(max(1e-19,minimum(jerror)), max(maximum(jerror),1e-4))
            title("Difference")
            ylabel(L"\Delta")
            axes_labels_opposite(:y)
        end
        csubplot(m,n,ks[2,3],nox=true) do
            j′error = abs.(j′-j′ref)
            loglog(x, j′error, rasterized=true)
            ylim(max(1e-19,minimum(j′error)), max(maximum(j′error),1e-4))
            ylabel(L"\Delta")
            axes_labels_opposite(:y)
        end
        csubplot(m,n,ks[3,3],nox=true) do
            yerror = abs.(y-yref)
            semilogx(x, log10.(yerror), rasterized=true)
            ylabel(L"\log_{10}(\Delta)")
            axes_labels_opposite(:y)
        end
        csubplot(m,n,ks[4,3]) do
            semilogx(x, log10.(abs.(y′-y′ref)), rasterized=true)
            xlabel(L"x")
            ylabel(L"\log_{10}(\Delta)")
            axes_labels_opposite(:y)
        end
    end
    savefig(joinpath(FIGDIR, "accuracy.svg"))
end

function continued_fractions()
    x = 10.0 .^ range(-2, stop=4, length=1_000)
    η = 10.0 .^ (-1:3)

    cfb = zeros(Float64, length(x))
    cfc1₊ = zeros(Float64, length(x), length(η))
    cfc2₊ = zeros(ComplexF64, length(x), length(η))
    cfc1₋ = zeros(Float64, length(x), length(η))
    cfc2₋ = zeros(ComplexF64, length(x), length(η))

    itersb = zeros(Int, length(x))
    itersc1₊ = zeros(Int, length(x), length(η))
    itersc2₊ = zeros(Int, length(x), length(η))
    itersc1₋ = zeros(Int, length(x), length(η))
    itersc2₋ = zeros(Int, length(x), length(η))

    λ = 0

    nx = length(x)
    nη = length(η)
    for (i,x) in enumerate(x)
        cfb[i],itersb[i] = bessel_fraction(x, λ)[1:2]
        for (j,η) in enumerate(η)
            cfc1₊[i,j],itersc1₊[i,j] = coulomb_fraction1(x, η, λ)[1:2]
            cfc2₊[i,j],itersc2₊[i,j] = coulomb_fraction2(x, η, λ, 1)[1:2]
            cfc1₋[i,j],itersc1₋[i,j] = coulomb_fraction1(x, -η, λ)[1:2]
            cfc2₋[i,j],itersc2₋[i,j] = coulomb_fraction2(x, -η, λ, 1)[1:2]
        end
    end

    cfigure("continued fractions", figsize=(9,6)) do
        csubplot(121) do
            loglog(x, itersc1₊)
            loglog(x, itersc2₊, "--")
            loglog(x, itersb, "k")
            text(x[1], 0.7itersc1₊[1,1], L"\eta=10^{-1}")
            text(x[1], 1.3itersc1₊[1,end], L"\eta=10^{3}", rotation=24)
            text(x[1], 0.25itersc2₊[1,1], L"\eta=10^{-1}", rotation=-65, horizontalalignment="left", verticalalignment="center")
            text(x[1], 0.8itersc2₊[1,end], L"\eta=10^{3}", rotation=-50, horizontalalignment="left", verticalalignment="center")
            title(L"\eta > 0")
            xlabel(L"x")
        end
        csubplot(122) do
            loglog(x, itersc1₋)
            loglog(x, itersc2₋, "--")
            loglog(x, itersb, "k")
            text(x[1], 0.7itersc1₋[1,1], L"\eta=-10^{-1}")
            text(x[1], 1.3itersc1₋[1,end], L"\eta=-10^{3}", rotation=33)
            text(x[1], 0.2itersc2₋[1,1], L"\eta=-10^{-1}", rotation=-65, horizontalalignment="left", verticalalignment="center")
            text(x[1], 0.7itersc2₋[1,end], L"\eta=-10^{3}", rotation=-55, horizontalalignment="left", verticalalignment="center")
            xlabel(L"x")
            title(L"\eta < 0")
            axes_labels_opposite(:y)
        end
    end
    savefig(joinpath(FIGDIR, "continued-fractions.svg"))
end

macro echo(expr)
    println(expr)
    :(@time $expr)
end

@info "Documentation plots"
mkpath(FIGDIR)
@echo simple_example()
@echo simple_coulomb_example()
@echo accuracy()
@echo continued_fractions()
