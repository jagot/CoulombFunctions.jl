using PyPlot
using Jagot.plotting
plot_style("ggplot")

using PyCall
Cycler = pyimport("cycler")
plt.rc("axes", prop_cycle=plt.rcParams["axes.prop_cycle"])

using SphericalBesselFunctions

include(joinpath(@__DIR__, "../test/reference.jl"))

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
    savefig("docs/src/figures/simple-example.svg")
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
            title("SphericalBesselFunctions.jl")
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
    savefig("docs/src/figures/accuracy.svg")
end

macro echo(expr)
    println(expr)
    :(@time $expr)
end

@info "Documentation plots"
mkpath("docs/src/figures")
@echo simple_example()
@echo accuracy()
