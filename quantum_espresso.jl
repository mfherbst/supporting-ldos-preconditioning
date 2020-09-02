using PyCall
using DFTK

# Setup ASE quantum espresso calculator
ENV["ASE_ESPRESSO_COMMAND"] = "/usr/bin/pw.x -in PREFIX.pwi > PREFIX.pwo"
ENV["ESPRESSO_PSEUDO"] = joinpath(@__DIR__, "qe_pseudos")

const PSEUDOMAP = Dict(
   "Si" => "Si.pbe-hgh.UPF",
   "O"  => "O.pbe-hgh.UPF",
   "Al" => "Al.pbe-hgh.UPF",
   "Cu" => "Cu.pbe-d-hgh.UPF",
   "Li" => "Li.pz-s-hgh.UPF",
   "Ga" => "Ga.pbe-d-hgh.UPF",
   "As" => "As.pbe-hgh.UPF",
   "H" => "H.pbe-hgh.UPF",
)

"""Ecut is in Hartree, tol is energy tolerance in Hartree, kspacing is in inverse Bohrs"""
function run_qe(lattice, atoms; Ecut=20, tol=1e-11, mixing_mode="local-TF",
                fileprefix="espresso", kspacing=0.3, smearing="gaussian",
                temperature=0, mixing_beta=0.7, mixing_ndim=10,
                maxiter=500, kwargs...)
    espresso = pyimport("ase.calculators.espresso")

    kgrid = kgrid_size_from_minimal_spacing(lattice, kspacing)
    @info "Autoselected kgrid_size", kgrid
    system = ase_atoms(lattice, atoms)
    calc = espresso.Espresso(
        label=fileprefix,
        input_dft="pbe",
        pseudopotentials=PSEUDOMAP,
        kpts=kgrid,
        ecutwfc=Ecut / DFTK.units.Ry,  # QE uses Rydbergs
        tstress=true,                  # Compute stresses
        tprnfor=true,                  # Compute forces
        mixing_mode=mixing_mode,
        mixing_beta=mixing_beta,
        conv_thr=tol / DFTK.units.Ry,
        occupations=ifelse(temperature > 0, "smearing", "fixed"),
        smearing=smearing,
        degauss=temperature / DFTK.units.Ry,
        electron_maxstep=maxiter,
        mixing_ndim=mixing_ndim,
        # verbosity="high",
        kwargs...
    )
    system.calc = calc

    system.get_potential_energy() * DFTK.units.eV
end

function parse_pwo(fileprefix)
    fn = isfile(fileprefix) ? fileprefix : fileprefix * ".pwo"
    iterations = []
    current_iter = nothing

    for line in readlines(fn)
        miter = match(r"^ *iteration # *([0-9]+)", line)
        if miter !== nothing  # New iteration starts
            if !isnothing(current_iter)
                push!(iterations, (; current_iter...))
            end
            current_iter = Dict{Symbol, Any}(:n_iter => parse(Int, miter[1]))
        end

        menergy = match(r"^ *total energy *= *([0-9.eE+-]+) *Ry", line)
        if menergy !== nothing
            current_iter[:total_energy] = DFTK.units.Ry * parse(Float64, menergy[1])
        end

        maccuracy = match(r"^ *estimated scf accuracy *< *([0-9.eE+-]+) *Ry", line)
        if maccuracy !== nothing
            current_iter[:accuracy] = DFTK.units.Ry * parse(Float64, maccuracy[1])
        end
    end

    iterations
end
