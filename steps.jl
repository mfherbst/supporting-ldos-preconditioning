include("builder.jl")
include("dielectric.jl")
include("dftk.jl")
include("quantum_espresso.jl")
using FFTW
using LinearAlgebra
using DataFrames
import JSON
import TimerOutputs


function setup_threading()
    FFTW.set_num_threads(4)
    BLAS.set_num_threads(4)
end


function random_logfile(prefix)
    logfile = "$(prefix)_$(@sprintf "%05d" rand(UInt16)).log"
    while isfile(logfile)
        logfile = "$(prefix)_$(@sprintf "%05d" rand(UInt16)).log"
    end
    logfile
end


function run_step_on_cases(payload, cases)
    spay = string(payload)
    header_printed = false

    # TODO Stopfile to stop iterations

    anymissing = false
    for dir in sort(collect(keys(cases)))
        config = cases[dir]

        if isfile(joinpath(config.directory, string(payload) * ".done"))
            continue
        end

        for i in 1:10  # To avoid infinite loop
            if !header_printed
                println("#####" * "#"^length(spay) * "#####")
                println("#--  $(string(payload))  --#")
                println("#####" * "#"^length(spay) * "#####")
                header_printed = true
            end

            logfile = random_logfile(joinpath(config.directory, string(payload)))
            println(basename(dir), "  ->  ", basename(dir) * "/" * basename(logfile))
            res = open(logfile, "w") do fp
                redirect_stdout(fp) do
                    payload(config)
                end
            end
            println("    ... $res\n")
            if res == :exhausted
                rm(logfile)
                anymissing || touch(joinpath(config.directory, string(payload) * ".done"))
                break
            elseif res == :missing
                anymissing = true
                rm(logfile)
            end
        end
    end
end


function firstrun(config)
    lockfile = joinpath(config.directory, "firstrun.running")
    firstprefix = joinpath(config.directory, "firstrun")
    (isfile(lockfile) || has_results(firstprefix)) && return :exhausted
    write(lockfile, read(`hostname`, String))

    lattice, atoms = load_system(config.case, config.n_repeat)
    firstα = get(config, :firstα, 0.5)
    localiser = config.construct_localiser(lattice, atoms)
    firstmixing = config.mixings[config.firstmixkey](firstα, localiser)

    try
        println("#\n#-- Firstrun $(config.case) n_repeat=$(config.n_repeat) -- " *
                "mixing=$firstmixing\n#")
        TimerOutputs.reset_timer!(DFTK.timer)
        run_dftk(lattice, atoms; fileprefix=firstprefix, mixing=firstmixing,
                 maxiter=100, config.kwargs_scf...)
        println(DFTK.timer)
    catch
        rm(joinpath(config.directory, "firstrun.done"), force=true)
        rethrow()
    finally
        rm(lockfile)
    end

    :done
end


function optimal_damping(config)
    lattice, atoms = load_system(config.case, config.n_repeat)
    localiser = config.construct_localiser(lattice, atoms)

    # Load firstrun results
    firstprefix = joinpath(config.directory, "firstrun")
    has_results(firstprefix) || return :missing
    firstmixing = config.mixings[config.firstmixkey](0.5, localiser)
    scfres = run_dftk(lattice, atoms; fileprefix=firstprefix, mixing=firstmixing,
                      config.kwargs_scf...)

    mixkey = nothing
    lockfile = nothing
    resultfile = nothing
    for mix in sort(collect(keys(config.mixings)))
        lockfile = joinpath(config.directory, mix * ".running")
        resultfile = joinpath(config.directory, mix * ".optimal_damping.json")
        (isfile(lockfile) || isfile(resultfile)) && continue
        write(lockfile, read(`hostname`, String))
        mixkey = mix
        break
    end
    mixkey === nothing && return :exhausted

    try
        mixing = config.mixings[mixkey](1.0, localiser)
        println("#\n#-- optimal_damping $(config.case) n_repeat=$(config.n_repeat) " *
                "-- mixing=$mixing\n#")
        fileprefix = joinpath(config.directory, "$(mixkey)_optimal_damping")

        TimerOutputs.reset_timer!(DFTK.timer)
        arnoldires = run_arnoldi_dielectric(scfres, mixing=mixing, force_lda=true,
                                            fileprefix=fileprefix, iterations=60)
        println(DFTK.timer)

        # ε = 1 - χ0 (vc + fxc)
        # In linear response regime: F(ρₙ) = (1 - ε) ρₙ
        # SCF Damped fixed-point scheme is
        #     ρₙ₊₁ = ρₙ + α (F(ρₙ) - ρₙ)
        #          = ρₙ + α (1 - ε - 1)ρₙ
        #          = (1 - αε) ρₙ
        # which converges iff
        #    σ(1 - αε) ⊂ (-1, 1)
        # but where things converge better the further the eigenvalues are from 1 and -1.
        # Optimising properly the problem
        #    min_α max(|1 - αε|, |1 + αε|))
        # shows the best α to give
        #    min σ(1 - αε) = max σ(1 - αε) ∈ (-1, 1)
        # or
        #    α_opt = 2 / (λ_max(ε) + λ_min(ε))

        eigenvalues  = real(arnoldires.λ)
        best_damping = 2 / (maximum(abs, eigenvalues) + minimum(abs, eigenvalues))

        # Since we are doing LDA here but our calculation is PBE there is really no
        # point in keeping more than 3 significant digits
        best_damping = round(best_damping, sigdigits=3)

        open(resultfile, "w") do fp
            save_dict = Dict(
                "eigenvalues" => real(arnoldires.λ),
                "residual_norms" => real(arnoldires.residual_norms),
                "condition_number" => maximum(abs, eigenvalues) / minimum(abs, eigenvalues),
                "optimal_damping" => best_damping,
            )
            JSON.print(fp, save_dict)
        end
    catch
        rm(joinpath(config.directory, "optimal_damping.done"), force=true)
        rethrow()
    finally
        rm(lockfile)
    end

    :done
end


function scf(config)
    lattice, atoms = load_system(config.case, config.n_repeat)
    localiser = config.construct_localiser(lattice, atoms)

    mixkey = nothing
    lockfile = nothing
    resultfile = nothing
    fileprefix = nothing
    for mix in sort(collect(keys(config.mixings)))
        lockfile = joinpath(config.directory, mix * ".running")
        fileprefix = joinpath(config.directory, mix)
        (isfile(lockfile) || has_results(fileprefix))  && continue
        write(lockfile, read(`hostname`, String))
        mixkey = mix
        break
    end
    mixkey === nothing && return :exhausted

    try
        # Load optimal damping results
        resultfile = joinpath(config.directory, mixkey * ".optimal_damping.json")
        isfile(resultfile) || return :missing
        αopt = open(JSON.parse, resultfile, "r")["optimal_damping"]

        mixing = config.mixings[mixkey](αopt, localiser)
        println("#\n#-- SCF $(config.case) n_repeat=$(config.n_repeat) -- " *
                "mixing=$mixing\n#")

        TimerOutputs.reset_timer!(DFTK.timer)
        run_dftk(lattice, atoms; fileprefix=fileprefix, mixing=mixing, config.kwargs_scf...)
        println(DFTK.timer)
    catch
        rm(joinpath(config.directory, "scf.done"), force=true)
        rethrow()
    finally
        rm(lockfile)
    end
    plot_convergence_summary(config)

    :done
end


function collect_minimal_total_energy(config::NamedTuple)
    collect_minimal_total_energy(config.directory, config.mixings)
end
function collect_minimal_total_energy(directory::AbstractString, mixings::Dict)
    # First pass: For DFTK runs collect mimimal total energy
    #       (The implementation in QE is slightly different
    #        resulting in slightly deviating energies)
    mintotal = 0
    for mixkey in sort(collect(keys(mixings)))
        resultfile = joinpath(directory, mixkey * ".json")
        isfile(resultfile) || continue

        energies = open(JSON.parse, resultfile, "r")["energies"]
        energies === nothing && continue
        isempty(energies) && continue

        mintotal = min(mintotal, minimum(energies[max(end-4, 1):end]))
    end

    mintotal
end


function plot_convergence_summary(config;
                                  iteration=30, tol=1e-10,
                                  outputfile=joinpath(config.directory, "summary.svg"))
    style = Dict(
        "None"          => (color=4, linestyle=:dot,   mark=:cross),
        "TF"            => (color=5, linestyle=:dot,   mark=:xcross),
        "Kerker"        => (color=2, linestyle=:solid, mark=:xcross),
        "RestaSiO2"     => (color=0, linestyle=:solid, mark=:xcross),
        "RestaGaAs"     => (color=1, linestyle=:dash,  mark=:xcross),
        "HybridSiO2"    => (color=0, linestyle=:solid, mark=:dtriangle),
        "HybridGaAs"    => (color=1, linestyle=:dash,  mark=:star4),
        "HybridVac"     => (color=4, linestyle=:dot,   mark=:star6),
        "Custom"        => (color=5, linestyle=:dash,  mark=:star8),
    )
    errors = Dict{String, Vector{Float64}}()

    # First pass: For DFTK runs collect mimimal total energy
    #       (The implementation in QE is slightly different
    #        resulting in slightly deviating energies)
    mintotal = collect_minimal_total_energy(config)

    # Add DFTK errors
    for mixkey in sort(collect(keys(config.mixings)))
        resultfile = joinpath(config.directory, mixkey * ".json")
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
    qe_fn = joinpath(config.directory, "TF/espresso.pwo")
    if isfile(qe_fn)
        errors["TF"] = [p.accuracy for p in parse_pwo(qe_fn)]
    end

    # If no data return
    isempty(errors) && return plot()

    ymax = 0
    maxiter = 0
    p = plot(legend=:topright, title=config.case)
    for mixkey in sort(collect(keys(errors)))
        dampfile = joinpath(config.directory, mixkey * ".optimal_damping.json")
        αopt = isfile(dampfile) ? open(JSON.parse, dampfile, "r")["optimal_damping"] : NaN
        mixkey == "TF" && (αopt = 0.7)  # QE calculations all use α=0.7

        maxiter = max(maxiter, length(errors[mixkey]))
        ymax = max(ymax, maximum(errors[mixkey]))
        plot!(p, errors[mixkey]; label=mixkey * " (α=$αopt)", style[mixkey]...)
    end
    xlims!(p, (0, min(iteration, maxiter)))
    ylims!(p, (tol, 2 * ymax))
    xaxis!(p, "Iteration")
    yaxis!(p, "Total energy absolute error", :log10)

    if isnothing(outputfile)
        p
    else
        savefig(p, outputfile)
    end
end


function dump_ldos(config, temperatures=[1e-3, 2e-3, 4e-3, 8e-3, 16e-3, 32e-3, 64e-3, 128e-3])
    lattice, atoms = load_system(config.case, config.n_repeat)
    localiser = config.construct_localiser(lattice, atoms)

    # Load firstrun results
    firstprefix = joinpath(config.directory, "firstrun")
    has_results(firstprefix) || return nothing
    firstmixing = config.mixings[config.firstmixkey](0.5, localiser)
    scfres = run_dftk(lattice, atoms; fileprefix=firstprefix, mixing=firstmixing,
                      config.kwargs_scf...)

    isnothing(temperatures) && (temperatures = [scfres.ham.basis.model.temperature])
    h5open(joinpath(config.directory, "ldos.hdf5"), "w") do h5f
        h5f["temperatures"] = temperatures
        h5f["temperature_scf"] = scfres.ham.basis.model.temperature
        h5f["ρ_real", "shuffle", (), "compress", 3] = real(scfres.ρ.real)
        for temp in temperatures
            stemp = @sprintf "%.3f" temp
            ldos = LDOS(scfres.εF, scfres.ham.basis, scfres.eigenvalues, scfres.ψ,
                        temperature=temp)
            h5f["ldos_$stemp", "shuffle", (), "compress", 3] = ldos
        end
    end
end


function collect_results(config; tol=1.6e-10)
    df = DataFrame(case=String[],
                   repeat=Int[],
                   size=Float64[],
                   mixing=String[],
                   damping=Union{Float64,Missing}[],
                   iterations=Union{Missing,Int}[],
                   condition_number=Union{Float64,Missing}[],
                   λmin=Union{Float64,Missing}[],
                   λmax=Union{Float64,Missing}[],
                   resmin=Union{Float64,Missing}[],
                   resmax=Union{Float64,Missing}[],
                  )

    for cfgkey in sort(collect(keys(config)))
        cfg        = config[cfgkey]
        minenergy  = collect_minimal_total_energy(cfg)
        lattice, _ = load_system(cfg.case, cfg.n_repeat)

        # Some made-up size function for now:
        size = norm(lattice * [0.0, 0.0, 1.0])

        for mixkey in sort(collect(keys(cfg.mixings)))
            scfdata  = Dict{String, Any}()
            dampdata = Dict{String, Any}()
            dampfile = joinpath(cfg.directory, mixkey * ".optimal_damping.json")
            scffile  = joinpath(cfg.directory, mixkey * ".json")

            isfile(dampfile) && (dampdata = open(JSON.parse, dampfile, "r"))
            isfile(scffile) && (scfdata  = open(JSON.parse, scffile, "r"))

            iterations = missing
            if "energies" in keys(scfdata)
                itfirst = findfirst(abs.(scfdata["energies"] .- minenergy) .≤ tol)
                iterations = isnothing(itfirst) ? typemax(Int) : itfirst
            end
            push!(df, (cfg.case,
                       cfg.n_repeat,
                       size,
                       mixkey,
                       get(dampdata, "optimal_damping", missing),
                       iterations,
                       get(dampdata, "condition_number", missing),
                       get(dampdata, "eigenvalues", [missing])[1],
                       get(dampdata, "eigenvalues", [missing])[end],
                       get(dampdata, "residual_norms", [missing])[1],
                       get(dampdata, "residual_norms", [missing])[end],
                 ))
        end
    end

    df
end
