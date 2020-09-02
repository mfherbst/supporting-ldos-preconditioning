using Plots
using Roots

function plot_dielectric_functions()
    χ0_metal(q, ktf=1) = -ktf^2 / (4π)
    χ0_dielectric(q, εr, C₀=1-εr, ktf=1) = C₀*q^2 / (4π * (1 - C₀*q^2/ktf^2))
    ε(χ0, q) = (1 - (4π)/q^2 * χ0(q))

    kerker(q, ktf=1) = q^2 / (ktf^2 + q^2)

    xs = range(1e-5, 1.5, length=10000)
    p = plot(xlims=(0, Inf), )
    plot!(p, x -> ε(χ0_metal, x),                    xs, label="Al",   ls=:solid, lw=2)
    plot!(p, x -> ε(q -> χ0_dielectric(q, 14.9), x), xs, label="GaAs", ls=:dash, lw=2)
    plot!(p, x -> ε(q -> χ0_dielectric(q,  2.4), x), xs, label="SiO₂", ls=:dashdot, lw=2)
    ylims!(p, (0, 16))
    xaxis!(p, "q")
    yaxis!(p, "ε(q)")

    p
end


function intersection_kerker_preconditioned()
    # Find intersection between Kerker-preconditioned and not preconditioned

    function ratio(qmin; qmax=100, prec=q -> 1)
        a = prec(qmin) * ε(q -> χ0_dielectric(q, 14.9), qmin)
        b = prec(qmax) * ε(q -> χ0_dielectric(q, 14.9), qmax)
        max(a, b) / min(a, b)
    end
    p = plot(xlims=(0, Inf), )
    plot!(p, q -> ratio(q),              xs, ls=:solid, label="None")
    plot!(p, q -> ratio(q, prec=kerker), xs, ls=:dash,  label="Kerker")
    ylims!(p, (0, 16))
    xaxis!(p, "q")
    yaxis!(p, "κ(q)")
    # savefig(p, "prec_dielectric_functions.pdf")

    qcut = Roots.find_zero(x -> ratio(x) - ratio(x, prec=kerker), [0.05, 0.1])
    println("qcut = $qcut")
    println("L    = $(2π / qcut)")
end
