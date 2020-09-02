include("../quantum_espresso.jl")
include("../builder.jl")
config = Dict{String, Any}()
include("config.jl")

function main()
    cfg = config[@__DIR__]
    kwargs_scf = cfg.kwargs_scf
    mixing_beta = 0.7

    lattice, atoms = load_system(cfg.case, cfg.n_repeat)
    fileprefix = joinpath(cfg.directory, "TF/espresso")
    run_qe(lattice, atoms, fileprefix=fileprefix, mixing_beta=mixing_beta,
           temperature=kwargs_scf.temperature, kspacing=kwargs_scf.kspacing,
           smearing="gaussian", maxiter=50)
end
