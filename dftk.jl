using DFTK
using Printf
using HDF5
using Dates
import FFTW
import JSON
using LinearAlgebra
import DFTK: HybridMixing
import DFTK: ScfDiagtol, ScfDefaultCallback

function write_scfres(h5gr, scfres; store_wfn=false)
    for symbol in (:ρin, :ρ)
        hasproperty(scfres, symbol) || continue
        object = getproperty(scfres, symbol)
        h5gr["$(symbol)_real", "shuffle", (), "compress", 3] = object.real
        h5gr["$(symbol)_fourier", "shuffle", (), "compress", 3] = object.fourier
    end

    for symbol in (:εF, )
        hasproperty(scfres, symbol) || continue
        object = getproperty(scfres, symbol)
        object === nothing && continue
        h5gr[string(symbol)] = object
    end

    for symbol in (:eigenvalues, :occupation)
        hasproperty(scfres, symbol) || continue
        h5gr[string(symbol)] = hcat(getproperty(scfres, symbol)...)
    end

    for key in keys(scfres.energies)
        h5gr["energy_$(key)"] = scfres.energies[key]
    end

    if store_wfn
        for symbol in (:ψ, )
            hasproperty(scfres, symbol) || continue
            object = getproperty(scfres, symbol)
            object === nothing && continue

            for ik in 1:length(object)
                h5gr["ψ_$(ik)_fourier", "shuffle", (), "compress", 3] = object[ik]
            end
        end
    end
end


function load_scfres(basis, fileprefix)
    has_results(fileprefix) || error("Need result in $fileprefix")

    h5open(fileprefix * ".hdf5", "r") do h5f
        res = h5f["result"]
        n_kpoints = length(basis.kpoints)
        n_bands = size(res["ψ_1_fourier"], 2)

        # Checks for consistency
        @assert size(res["ρ_real"]) == basis.fft_size
        @assert size(res["occupation"]) == (n_bands, n_kpoints)
        @assert size(res["eigenvalues"]) == (n_bands, n_kpoints)

        # Read state
        εF = read(res, "εF")
        ρ = from_real(basis, res["ρ_real"][:, :, :])
        ψ = Matrix{ComplexF64}[]
        occupation = Vector{Float64}[]
        eigenvalues = Vector{Float64}[]
        for ik = 1:n_kpoints
            n_fft = length(G_vectors(basis.kpoints[ik]))
            @assert size(res["ψ_$(ik)_fourier"]) == (n_fft, n_bands)
            push!(ψ, res["ψ_$(ik)_fourier"][:, :])

            push!(occupation, res["occupation"][:, ik])
            push!(eigenvalues, res["eigenvalues"][:, ik])
        end

        # Construct Hamiltonian and compute energies
        energies, ham = energy_hamiltonian(basis, ψ, occupation;
                                           ρ=ρ, eigenvalues=eigenvalues, εF=εF)

        (ham=ham, basis=basis, energies=energies, converged=true, ρ=ρ,
         eigenvalues=eigenvalues, occupation=occupation, εF=εF, n_iter=nothing, ψ=ψ)
    end
end


function has_results(fileprefix)
    !isfile(fileprefix * ".hdf5") && return false
    try
        return h5open(fileprefix * ".hdf5", "r") do h5f
            exists(h5f, "result")
        end
    catch
    end
    false
end


function run_dftk(lattice, atoms; Ecut=20, tol=1e-11, mixing=DFTK.SimpleMixing(),
                  maxiter=50, kspacing=0.3, smearing=Smearing.None(),
                  temperature=0, fileprefix="dftk",
                  determine_bands=nelec -> nelec // 2 + sqrt(nelec // 2),
                  stepsize=10, randomise_guess=false, kgrid=nothing,
                  supersampling=1.5, ensure_odd_fft_size=false,
                  find_occupation=find_occupation, random_ampl=1e-2,
                  diagtol_max=1e-3,
                 )
    model = model_DFT(lattice, atoms, [:gga_x_pbe, :gga_c_pbe],
                      smearing=smearing, temperature=temperature)
    isnothing(kgrid) && (kgrid = kgrid_size_from_minimal_spacing(lattice, kspacing))

    fft_size = determine_grid_size(model, Ecut, supersampling=supersampling)
    ensure_odd_fft_size && (fft_size = Vec3(fft_size .+ iseven.(fft_size)))
    basis = PlaneWaveBasis(model, Ecut, kgrid=kgrid, kshift=[0, 0, 0], fft_size=fft_size)
    n_bands = ceil(Int, determine_bands(model.n_electrons))

    if has_results(fileprefix)
        @info "Loading existing results instead of SCF"
        return load_scfres(basis, fileprefix)
    end
    n_blas = Int(ccall((BLAS.@blasfunc(openblas_get_num_threads), BLAS.libblas), Cint, ()))

    println()
    println("hostname       = $(read(`hostname`, String))")
    println("started on     = $(Dates.now())")
    println("julia threads  = $(Threads.nthreads())")
    println("BLAS threads   = $(n_blas)")
    println()
    println("temperature    = $(basis.model.temperature)")
    println("smearing       = $(basis.model.smearing)")
    if isdiag(lattice)
        println("diag(lattice)  = $(round.(diag(lattice), sigdigits=4))")
    else
        println("lattice        = $(round.(lattice, sigdigits=4))")
    end
    println("Ecut           = $Ecut")
    println("fft_size       = $(basis.fft_size)")
    println("kgrid          = $kgrid")
    println("irreducible k  = $(length(basis.kpoints))")
    println("n_bands        = $(n_bands)")
    println("n_electr       = $(basis.model.n_electrons)")
    println("mixing         = $mixing")
    flush(stdout)

    # Collection of data to dump later
    save_dict = Dict(
        "temperature"  => basis.model.temperature,
        "smearing"     => string(basis.model.smearing),
        "diag_lattice" => diag(lattice),
        "lattice"      => lattice,
        "Ecut"         => Ecut,
        "fft_size"     => basis.fft_size,
        "kgrid"        => kgrid,
        "n_kirred"     => length(basis.kpoints),
        "n_bands"      => n_bands,
        "n_electr"     => basis.model.n_electrons,
        "mixing"       => string(mixing),
        "α"            => mixing.α,
        "energies"     => Float64[],
    )

    # Empty the existing HDF5 file
    isdir(dirname(fileprefix)) || mkdir(dirname(fileprefix))
    close(h5open(fileprefix * ".hdf5", "w"))

    default_callback = ScfDefaultCallback()
    function callback(info)
        default_callback(info)
        flush(stdout)
        key = info.stage == :finalize ? "result" : "neval_$(info.n_iter)"
        store_wfn = info.stage == :finalize
        h5open(fileprefix * ".hdf5", "r+") do h5f
            write_scfres(g_create(h5f, key), info, store_wfn=store_wfn)
        end
        push!(save_dict["energies"], info.energies.total)
    end

    ρ0 = guess_density(basis)
    if randomise_guess
        # Add randomness and adjust electron count
        ρ0.real .= abs.(ρ0.real .+ random_ampl * rand(size(ρ0.real)))
        ρ0.real .*= basis.model.n_electrons ./ (sum(ρ0.real) ./ length(ρ0.real) .* basis.model.unit_cell_volume)
    end
    scfres = self_consistent_field(basis, ρ=ρ0, tol=tol, mixing=mixing,
                                   callback=callback,
                                   determine_diagtol=ScfDiagtol(ratio_ρdiff=0.04,
                                                                diagtol_max=diagtol_max,
                                                                diagtol_min=1e-12),
                                   solver=scf_nlsolve_solver(stepsize),
                                   maxiter=maxiter, n_bands=n_bands,
                                   find_occupation=find_occupation,
                                   enforce_symmetry=true)

    println(scfres.energies)
    flush(stdout)

    open(fp -> JSON.print(fp, save_dict), fileprefix * ".json", "w")
    scfres
end
