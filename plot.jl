include("postprocessing.jl")
include("anderson.jl")
include("dielectric_functions.jl")
using Plots.PlotMeasures


function setup_plots()
    gr()
    default(size=tuple(Int.(ceil.(0.75 .* [600, 400]))...),
            guidefontsize=10, grid=false)
end


function dump_latex_tables()
    write(joinpath(@__DIR__, "table_bulk.tex"),  table_bulk(df,  latex=true))
    write(joinpath(@__DIR__, "table_mixed.tex"), table_mixed(df, latex=true))
end


function plot_modes(case="Al_repeat_40", colorbar=false)
    h5f = h5open(joinpath(@__DIR__, case, "None_optimal_damping.hdf5"))

    # Read modes and average over X axes
    xlen = size(h5f["arnoldi/X"], 1)
    modes = [
        sum(h5f["arnoldi/X"][:, :, :,   i], dims=1)[1, :, :] / xlen
        for i in 1:size(h5f["arnoldi/X"], 4)
    ]
    λs = h5f["arnoldi/λ"][:]

    data = [
            modes[end-0];
            fill(0, 10, size(modes[end], 2));
            modes[end-1];
            fill(0, 10, size(modes[end], 2));
            modes[end-2];
            fill(0, 10, size(modes[end], 2));
            modes[end-3];
            fill(0, 10, size(modes[end], 2));
            modes[end-4];
           ]

    p = heatmap(data, c=:RdBu_11, aspect_ratio=1, showaxis=false, ticks=false,
                margin=0mm, clim=(-0.006, 0.006), grid=false, legend=colorbar)
    yaxis!(p, raw"$y$ direction")
    xaxis!(p, raw"$z$ direction")

    # Draw black border
    w = size(modes[end], 1) + 0.5
    h = size(modes[end], 2) + 0.5
    xs = [0.5, 0.5, h,   h, 0.5]
    ys = [0.5,   w, w, 0.5, 0.5]
    for a in 0:4
        plot!(p, xs, ys .+ a .* (w .+ 9.5), color=:black, legend=false, lw=0.5)
    end

    # Determine atomic (y, z) coordinates
    system, _, repeat = split(case, "_")
    lattice, atoms = load_system(system, parse(Int, repeat))
    positions = vcat(last.(atoms)...)
    cleanup(x) = x > 0.98 ? x - 1 : x
    atom_yzs = [(cleanup(pos[3]) * size(modes[end], 2),
                 cleanup(pos[2]) * size(modes[end], 1)) for pos in positions]
    for a in 0:4
        scatter!(p, first.(atom_yzs) .+ 1,
                 last.(atom_yzs) .+ a .* (w .+ 9.5),
                 m=:circle, color=:black, markersize=1.4, markeralpha=0.5,
                 markerstrokealpha=0.5)
    end

    p
end


function plot_ldos_density_2d(case; clim=(0, 0.15), aspect_ratio=2)
    h5f = h5open(joinpath(@__DIR__, case, "ldos.hdf5"))
    stemp = @sprintf "%.3f" read(h5f, "temperature_scf")

    # Average over X axes
    xlen = size(h5f["ldos_$stemp"], 1)
    data = [
        sum(h5f["ldos_$stemp"][:, :, :], dims=1)[1, :, :] / xlen;
        fill(0, 10, size(h5f["ldos_$stemp"], 3));
        sum(h5f["ρ_real"][:, :, :], dims=1)[1, :, :] / xlen;
    ]

    p = heatmap(data, c=cgrad(:devon, rev=true), aspect_ratio=aspect_ratio,
                showaxis=false, ticks=false, size=(450, 200),
                margin=0mm, grid=false, clim=clim)
    yaxis!(p, raw"$y$ direction")
    xaxis!(p, raw"$z$ direction")

    # Draw black border
    w = size(h5f["ldos_$stemp"], 2) + 0.5
    h = size(h5f["ldos_$stemp"], 3) + 0.5
    xs = [0.5, 0.5, h,   h, 0.5]
    ys = [0.5,   w, w, 0.5, 0.5]
    plot!(p, xs, ys, color=:black, legend=false, lw=0.5)
    plot!(p, xs, ys .+ w .+ 9.5,
          color=:black, legend=false, lw=0.5)

    p
end


function plot_ldos_density_1d(case; boxes=[], temperatures=nothing, flip=false,
                              legend=:topright, letter=nothing)
    h5f = h5open(joinpath(@__DIR__, case, "ldos.hdf5"))
    isnothing(temperatures) && (temperatures = [read(h5f, "temperature_scf")])

    # Average over X axes
    xlen = size(h5f["ρ_real"], 1)
    ylen = size(h5f["ρ_real"], 2)
    rho_avg  = sum(sum(h5f["ρ_real"][:, :, :], dims=1), dims=2)[1, 1, :] / xlen / ylen;

    bmin = Inf
    bmax = -Inf
    ldos_avg = Vector{Any}(undef, length(temperatures))
    for (i, temp) in enumerate(temperatures)
        stemp = @sprintf "%.3f" temp
        ldos_avg[i] = sum(sum(h5f["ldos_$stemp"][:, :, :],
                              dims=1), dims=2)[1, 1, :] / xlen / ylen;
        bmin = min(bmin, min(minimum(ldos_avg[i]), minimum(rho_avg)))
        bmax = max(bmax, max(maximum(ldos_avg[i]), maximum(rho_avg)))
    end

    # Determine atomic positions for tick placement
    system, _, repeat = split(case, "_")
    lattice, atoms = load_system(system, parse(Int, repeat))
    positions = vcat(last.(atoms)...)
    ticks = [pos[3] * size(h5f["ρ_real"], 3) for pos in positions]

    p = plot(;legend=legend, xticks=ticks, xformatter=x -> "", grid=false, )
    for (s, e) in boxes
        c=:grey85
        plot!(p,
              [s,    s,    e,    e,    s] .* length(rho_avg),
              [bmin, bmax, bmax, bmin, bmin],
              fill=(0, c), label="", line=c)
    end
    lss = append!([:solid, :dot, :dash, :dashdot], fill(:solid, length(temperatures)))
    cols = [3, :orangered1, 4, 5, 6, 7, 8, 9]
    for (i, temp) in reverse(collect(enumerate(temperatures)))
        stemp = @sprintf "%.3f" temp
        label = "LDOS" * (length(temperatures) > 1 ? " (T=$stemp)" : "")
        plot!(p, ldos_avg[i], label=label, color=cols[i], ls=lss[i], lw=1.8)
    end
    plot!(p, rho_avg, label="ρ", color=1, ls=:solid, lw=3.5)
    xaxis!(p, "z", xflip=flip)

    if !isnothing(letter)
        start = flip ? length(rho_avg) - 5 : 3
        annotate!(p, start, 1.05bmax, Plots.text(letter, 10, :left))
    end

    p
end

const NAMEMAP = Dict(
    "RestaSiO2"  => raw"Dielectric (SiO₂)",
    "RestaGaAs"  => raw"Dielectric (GaAs)",
    "HybridGaAs" => "LDOS+Dielectric",
    "HybridVac"  => "LDOS",
    "Custom"     => "Localised",
    "TF"         => "TFW",
)

function plot_convergence(case; namemap=NAMEMAP, iteration=50, tol=1e-10,
                          preconditioners=["None", "Kerker"],
                          ylabel="Total energy absolute error",
                          show_alpha=false, ymax=1e3, letter=nothing)
    style = Dict(
        # both
        "None"          => (color=4, linestyle=:dot,   mark=:cross, markerstrokecolor=4),
        "Kerker"        => (color=2, linestyle=:solid, mark=:xcross, markerstrokecolor=2),
        # bulk
        "RestaSiO2"     => (color=0, linestyle=:solid, mark=:utriangle, markerstrokecolor=0),
        "RestaGaAs"     => (color=1, linestyle=:dash,  mark=:dtriangle, markerstrokecolor=1),
        # compare
        "HybridVac"     => (color=1, linestyle=:solid, mark=:utriangle, markerstrokecolor=1),
        "TF"            => (color=0, linestyle=:dash,  mark=:dtriangle, markerstrokecolor=0),
    )

    # First pass: For DFTK runs collect mimimal total energy
    directory = joinpath(@__DIR__, case)
    mintotal = collect_minimal_total_energy(directory, style)

    # Add DFTK errors
    errors = Dict{String, Vector{Float64}}()
    for mixkey in preconditioners
        resultfile = joinpath(directory, mixkey * ".json")
        isfile(resultfile) || continue

        energies = open(JSON.parse, resultfile, "r")["energies"]
        energies === nothing && continue
        values = abs.(energies .- mintotal)
        idx = findfirst(values .< 1e-16)
        isnothing(idx) || (values = values[1:idx-1])
        isempty(values) && continue

        errors[mixkey] = values
    end

    # Add Quantum-Espresso errors
    qe_fn = joinpath(directory, "TF/espresso.pwo")
    if isfile(qe_fn) && "TF" in preconditioners
        errors["TF"] = [p.accuracy for p in parse_pwo(qe_fn)]
    end

    # If no data return
    isempty(errors) && return plot()

    maxiter = 0
    p = plot(legend=:topright, grid=false)
    for mixkey in preconditioners
        dampfile = joinpath(directory, mixkey * ".optimal_damping.json")
        αopt = isfile(dampfile) ? open(JSON.parse, dampfile, "r")["optimal_damping"] : NaN
        mixkey == "TF" && (αopt = 0.7)  # QE calculations all use α=0.7

        maxiter = max(maxiter, length(errors[mixkey]))
        label = get(NAMEMAP, mixkey, mixkey)
        show_alpha && (label=label * " (α=$αopt)")
        plot!(p, errors[mixkey]; label=label, markersize=4, style[mixkey]..., lw=1.5)
    end
    xlims!(p, (0, min(iteration, maxiter)))
    ylims!(p, (tol, ymax))
    yaxis!(p, ylabel, :log10)
    xaxis!(p, "Iteration")
    if !isnothing(letter)
        pos = ymax / 10 * 9
        plot!(p, [1.1, 12], [pos, pos], label="", line=(0.3, :white), lw=22)
        annotate!(p, 1.2, pos, Plots.text(letter, 10, :left))
    end

    p
end


function main()
    setup_plots()

    for (letter, sys) in [("(b) Al+vacuum", "AlVac_repeat_10"), ("(a) Al", "Al_repeat_12"),
                          ("(c) GaAs", "GaAs_repeat_20")]
        savefig(plot_anderson("$(sys)_anderson", letter=letter), "$(sys)_anderson.pdf")
    end

    # Plots of dielectric functions
    savefig(plot_dielectric_functions(), "dielectric_functions.pdf")
    pgfplotsx(); savefig(plot_modes(), "aluminium_modes.pdf"); gr()

    # LDOS / density plots
    # savefig(plot_ldos_density_2d("AlSiO2H_repeat_20"),              "AlSiO2_ldos.pdf")
    # savefig(plot_ldos_density_2d("AlVac_repeat_20", clim=(0, 0.2)), "AlVac_ldos.pdf")
    savefig(plot_ldos_density_1d("AlSiO2H_repeat_20", boxes=[(0.395, 0.995),],
                                 temperatures=[0.001, 0.032], flip=true,
                                 letter="(b) Al+SiO₂"),
            "AlSiO2_ldos_1d.pdf")
    savefig(plot_ldos_density_1d("AlVac_repeat_20", boxes=[(0, 0.516), (0.99, 1.0)],
                                 temperatures=[0.001, 0.032], letter="(a) Al+vacuum"),
            "AlVac_ldos_1d.pdf")

    # Convergence plots
    ldos_vs_tf_precs = ["None", "HybridVac", "Kerker", "TF"]
    savefig(plot_convergence("AlVac_repeat_20", preconditioners=ldos_vs_tf_precs,
                             ylabel="Estimated SCF error", letter="(a) Al+vacuum"),
            "AlVac_repeat_20.pdf")
    savefig(plot_convergence("AlSiO2H_repeat_20", preconditioners=ldos_vs_tf_precs,
                             ylabel="Estimated SCF error", letter="(b) Al+SiO₂"),
            "AlSiO2H_repeat_20.pdf")

    bulk_precs = ["None", "RestaSiO2", "RestaGaAs", "Kerker"]
    savefig(plot_convergence("SiO2_repeat_39", preconditioners=bulk_precs, letter="(a) SiO₂"),
            "SiO2_repeat_39.pdf")
    savefig(plot_convergence("GaAs_repeat_40", preconditioners=bulk_precs, letter="(b) GaAs"),
            "GaAs_repeat_40.pdf")
    savefig(plot_convergence("Al_repeat_40", preconditioners=bulk_precs, letter="(c) Al"),
            "Al_repeat_40.pdf")
end


(abspath(PROGRAM_FILE) == @__FILE__) && main()
