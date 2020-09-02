include("steps.jl")


config = Dict{String, Any}()
for d in readdir()
    if isdir(d) && isfile(joinpath(d, "config.jl"))
        include(joinpath(d, "config.jl"))
    end
end


function main()
    setup_threading()

    run_step_on_cases(firstrun, config)
    run_step_on_cases(optimal_damping, config)
    run_step_on_cases(scf, config)
end
