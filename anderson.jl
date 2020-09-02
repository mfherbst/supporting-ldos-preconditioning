using Plots
import JSON

function plot_anderson(folder; tol=1e-10, iteration=40, ymax=1e2, letter=nothing)
    mixing = "None"
    isfile(joinpath(folder, mixing * ".optimal_damping.json")) || (mixing = "Kerker")

    energies = Dict(m => Float64.(open(JSON.parse,
                                       joinpath(folder, "anderson_$m.json"))["energies"])
                    for m in (0, 5, 10))
    mintotal = minimum(minimum, values(energies))
    mintotal -= tol / 100  # To avoid zero error (and thus -Inf on a log scale)

    errors = Dict{Int, Vector{Float64}}()
    for m in sort(collect(keys(energies)))
        values = abs.(energies[m] .- mintotal)
        idx = findfirst(values .< 1e-16)
        isnothing(idx) || (values = values[1:idx-1])
        isempty(values) && continue
        errors[m] = values
    end

    maxiter = 0
    ms = [:x, :+, :utriangle, :dtriangle]
    p = plot(legend=:topright, grid=false)
    for (i, m) in enumerate(sort(collect(keys(errors))))
        maxiter = max(maxiter, length(errors[m]))
        plot!(p, errors[m]; label="Anderson m=$m", m=ms[i], markerstrokecolor=i,
              linewidth=2)
    end

    fn = joinpath(folder, "$mixing.optimal_damping.json")
    κ = open(JSON.parse, fn)["condition_number"]
    rate_simple = (κ - 1) / (κ + 1)                # Convergence rate for simple mixing
    rate_anderson = (sqrt(κ) - 1) / (sqrt(κ) + 1)  # Convergence rate for anderson

    # Square because energy converges quadratically (as opposed to density residual)
    x = collect(1:0.1:maxiter)
    plot!(p, x, errors[0][1]  * rate_simple.^2x, label="Simple rate", ls=:dash,
          color=:black, linewidth=2)
    plot!(p, x, errors[10][1] * rate_anderson.^2x, label="Anderson rate", ls=:dot,
          color=:black, linewidth=2)

    xlims!(p, (0, iteration))
    ylims!(p, (tol, ymax))
    xaxis!(p, "Iteration")
    yaxis!(p, "Total energy absolute error", :log10)
    if !isnothing(letter)
        plot!(p, [1.1, 8], [18, 18], label="", line=(0.5, :white), lw=15)
        annotate!(p, 1.2, 1.8e1, Plots.text(letter, 10, :left))
    end

    p
end
