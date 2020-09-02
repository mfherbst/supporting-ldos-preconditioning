using PyCall
using DFTK
import JSON
using LinearAlgebra


function ase_bulk(args...; kwargs...)
    ase_build = pyimport("ase.build")
    if first(args) == "SiO2"
        # Lattice parameters for α-quartz from
        # https://www.materialsproject.org/materials/mp-549166/
        # All lengths in Ångström (used by ASE)
        a = b = 4.98530291
        c = 5.47852971
        α = β = 90
        γ = 120

        silicon = [  # Fractional coordinates
            [0.000000  0.471138  0.833333],
            [0.471138  0.000000  0.166666],
            [0.528862  0.528862  0.500000],
        ]
        oxygen = [  # Fractional coordinates
            [0.148541  0.734505  0.046133],
            [0.265495  0.414036  0.712800],
            [0.585964  0.851459  0.379466],
            [0.734505  0.148541  0.953867],
            [0.414036  0.265495  0.287200],
            [0.851459  0.585964  0.620534],
        ]
        symbols = append!(fill("Si", 3), fill("O", 6))
        pyimport("ase").Atoms(symbols=symbols, cell=[a, b, c, α, β, γ], pbc=true,
                              scaled_positions=vcat(vcat(silicon, oxygen)...))
    elseif first(args) == "GaAs"
        # Lattice parameters for ZnS-type GaAs from
        # Haynes, William M., ed. (2011). CRC Handbook of Chemistry and Physics
        # (92nd ed.). CRC Press. ISBN 978-1439855119.
        a = 5.65315  # for zinc-blende cell
        ase_build.bulk("GaAs", "zincblende"; a=a, kwargs...)
    elseif first(args) == "AlN"
        a = b = 3.12858814
        c = 5.01695500
        α = β = 90
        γ = 120
        aluminium = [
            [1/3  2/3  0.99928700],
            [2/3  1/3  0.49928700],
        ]
        nitrogen = [
            [1/3  2/3  0.38071300],
            [2/3  1/3  0.88071300],
        ]
        symbols = ["Al", "Al", "N", "N"]
        pyimport("ase").Atoms(symbols=symbols, cell=[a, b, c, α, β, γ], pbc=true,
                              scaled_positions=vcat(vcat(aluminium, nitrogen)...))
    else
        ase_build.bulk(args...; kwargs...)
    end
end


function bulk(args...; kwargs...)
    atoms_ase = ase_bulk(args...; kwargs...)
    atoms = load_atoms(atoms_ase)
    atoms = [ElementPsp(el.symbol, psp=load_psp(el.symbol, functional="pbe")) => position
             for (el, position) in atoms]
    load_lattice(atoms_ase), atoms
end


function randomised_supercell(name; n_repeat=5)
    kwargs = ()
    if name == "Al"
        factor = 1/4  # 4 Al per repeat
        kwargs = (cubic=true, )
        amplitude = 1e-2
    elseif name == "GaAs"
        factor = 1    # One Ga and one As per repeat
        amplitude = 1e-2
    elseif name == "SiO2"
        factor = 1/3  # One Si and two O per repeat
        amplitude = 1e-2
    else
        error("Not implemented: $system")
    end
    @assert isinteger(n_repeat * factor)
    n_repeat = Int(n_repeat * factor)


    # Make supercell
    prim = ase_bulk(name; kwargs...)
    super = pyimport("ase.build").make_supercell(prim,
                                                 Array(Diagonal([1, 1, n_repeat])))

    # Introduce some random noise
    cartesian_pos = super.get_positions()
    ampl_Angstr = amplitude / DFTK.units.Å
    # noise = ampl_Angstr * randn(size(cartesian_pos))
    noise = 4ampl_Angstr .* (-0.5 .+ 0.5*rand(size(cartesian_pos)...))
    super.set_positions(cartesian_pos .+ noise)

    # Return DFTK datastructures
    atoms = [ElementPsp(el.symbol, psp=load_psp(el.symbol, functional="pbe")) => position
             for (el, position) in load_atoms(super)]
    load_lattice(super), atoms
end


function ase_surface(identifier; n_repeat=5, distance=nothing)
    # n_repeat gives the number of repeats in the surfaces
    # (i.e. n_repeat == 5 means 5 cells one material and 5 cells the other)
    ase_build = pyimport("ase.build")
    repeat = Array(Diagonal([1, 1, n_repeat]))  # repeat along z-dirn

    function get_distance(surface, user_set)  # in Å
        if isnothing(user_set)
            return maximum(abs, hcat(surface.cell...))
        else
            return user_set / DFTK.units.Å
        end
    end

    if identifier == "AlGaAs"
        gaas_110 = ase_lego("GaAs", (1, 1, 0), n_repeat=n_repeat)
        al_100 = ase_lego("Al", (1, 0, 0), n_repeat=n_repeat)
        return ase_build.stack(gaas_110, al_100 * (2, 1, 1), fix=0, maxstrain=1.7, distance=3.0)
    elseif identifier == "AlSiO2H"
        sio2h_110 = ase_lego("SiO2H", (1, 1, 0), cubic=false, n_repeat=n_repeat)
        ase_build.add_vacuum(sio2h_110, 1.0)
        sio2h_110.center()

        al_100 = ase_lego("Al", (1, 0, 0), n_repeat=n_repeat) * (2, 1, 1)
        ase_build.add_vacuum(al_100, 1.0)
        al_100.center()

        return ase_build.stack(sio2h_110, al_100, fix=0, maxstrain=1.6, distance=2.5)
    elseif identifier == "AlVac"
        al_100 = ase_lego("Al", (1, 0, 0), n_repeat=n_repeat)
        ase_build.add_vacuum(al_100, get_distance(al_100, distance))
        return al_100
    elseif identifier == "GaAsVac"
        gaas_110 = ase_lego("GaAs", (1, 1, 0), n_repeat=n_repeat)
        ase_build.add_vacuum(gaas_110, get_distance(gaas_110, distance))
        return gaas_110
    elseif identifier == "GaAsSiO2H"
        gaas_110 = ase_lego("GaAs", (1, 1, 0), n_repeat=n_repeat)
        sio2h_110 = ase_lego("SiO2H", (1, 1, 0), cubic=false, n_repeat=n_repeat)
        return ase_build.stack(gaas_110, sio2h_110, fix=0, maxstrain=0.7, distance=2.5)
    elseif identifier == "SiO2HVac"
        sio2h_110 = ase_lego("SiO2H", (1, 1, 0), cubic=false, n_repeat=n_repeat)
        ase_build.add_vacuum(sio2h_110, get_distance(sio2h_110, distance))
        return sio2h_110
    elseif identifier == "AlGaAsSiO2H"
        gaas_110  = ase_lego("GaAs", (1, 1, 0), n_repeat=n_repeat)
        sio2h_110 = ase_lego("SiO2H", (1, 1, 0), cubic=false, n_repeat=n_repeat)
        interm = ase_build.stack(gaas_110, sio2h_110, fix=0, maxstrain=0.7, distance=2.5)

        al_100 = ase_lego("Al", (1, 0, 0), n_repeat=n_repeat) * (2, 1, 1)
        return ase_build.stack(interm, al_100, fix=0, maxstrain=1.7, distance=2.5)
    else
        error("Unknown surface identifier: $identfier")
    end
end


function surface(args...; kwargs...)
    atoms_ase = ase_surface(args...; kwargs...)
    atoms = load_atoms(atoms_ase)
    atoms = [ElementPsp(el.symbol, psp=load_psp(el.symbol, functional="pbe")) => position
             for (el, position) in atoms]
    load_lattice(atoms_ase), atoms
end


function store_system(lattice, atoms, name, repeat)
    atnums = Int[]
    pspfiles = String[]
    positions = Vector{Float64}[]
    for (elem, pos) in atoms
        append!(atnums, fill(elem.Z, length(pos)))
        append!(pspfiles, fill(elem.psp.identifier, length(pos)))
        append!(positions, pos)
    end

    if !isdir(joinpath(@__DIR__, "systems"))
        mkdir(joinpath(@__DIR__, "systems"))
    end
    open(joinpath(@__DIR__, "systems/$(name)_$repeat.json"), "w") do fp
        save_dict = Dict(
            "lattice" => Array(lattice),
            "atnums" => atnums,
            "pspfiles" => pspfiles,
            "positions" => positions,
        )
        JSON.print(fp, save_dict)
    end
    nothing
end


function load_system(name, repeat)
    store = joinpath(@__DIR__, "systems/$(name)_$repeat.json")
    if !isfile(store)
        if name in ("Al", "SiO2", "GaAs")
            lattice, atoms = randomised_supercell(name, n_repeat=repeat)
        elseif name == "SiO2H"
            error("Sorry cannot build passivated SiO2")
        else
            lattice, atoms = surface(name, n_repeat=repeat)
        end
        store_system(lattice, atoms, name, repeat)
    end

    saved_dict = open(JSON.parse, store, "r")
    lattice = Mat3(hcat(saved_dict["lattice"]...))
    atoms = map(unique(saved_dict["atnums"])) do atnum
        pspid = saved_dict["pspfiles"][findfirst(atnum .== saved_dict["atnums"])]
        elem = ElementPsp(atnum, psp=load_psp(pspid))
        coords = Vec3.(saved_dict["positions"][atnum .== saved_dict["atnums"]])
        elem => coords
    end
    lattice, atoms
end


function ase_lego(system, miller; cubic=true, n_repeat=1)
    ase_build = pyimport("ase.build")
    ase_io = pyimport("ase.io")

    extra = ""
    if cubic
        extra = "_cubic"
    end
    mstr = prod(string, miller)

    store = joinpath(@__DIR__, "systems/ase_lego_$(system)$(extra)_$(mstr)_$(n_repeat).xyz")
    if !isfile(store)
        atoms = ase_bulk(system, cubic=cubic)

        factor = 1.0
        if system == "GaAs"
            factor = 0.5  # Two Ga and two As per layer
        elseif system == "SiO2"
            factor = 0.5  # Three Si and 6 O per layer
        elseif system == "Al"
            factor = 0.5  # 2 Al per layer
        else
            error("Not implemented: $system")
        end
        @assert isinteger(n_repeat * factor)
        n_repeat = Int(n_repeat * factor)
        atoms = ase_build.surface(atoms, miller, n_repeat, 0; periodic=true)
        ase_io.write(store, atoms)
    end

    ase_io.read(store)
end
