include("steps.jl")
using PrettyTables
import ColorSchemes

if !isdefined(Main, :config) || isnothing(config)
    config = Dict{String, Any}()
    for d in readdir()
        if isdir(d) && isfile(joinpath(d, "config.jl"))
            include(joinpath(d, "config.jl"))
        end
    end
end

df = collect_results(config)
case_cleanup = Dict(
    "SiO2HVac"    => "SiO2 / vac",
    "GaAsVac"     => "GaAs / vac",
    "AlVac"       => "Al / vac",
    "GaAsSiO2H"   => "GaAs / SiO2",
    "AlSiO2H"     => "Al / SiO2",
    "AlGaAs"      => "Al / GaAs",
    "AlGaAsSiO2H" => "Al / GaAs / SiO2",
)
mixing_cleanup = Dict(
    "Resta"      => "Dielectric",
    "RestaGaAs"  => "Dielectric (GaAs)",
    "RestaSiO2"  => "Dielectric (SiO2)",
    "HybridVac"  => "LDOS",
    "HybridGaAs" => "LDOS + Dielectric"
)


function ft_sigdigits(columns, sigdigits=3)
    function inner(v, i, j)
        ismissing(v) && return v
        if j in columns && applicable(round, v)
            return round(v, sigdigits=3)
        end
        v
    end
end
function ft_iterations(columns)
    function ft_iterations(v, i, j)
        ismissing(v) && return v
        if j in columns && v == typemax(Int)
            return ">50"
        end
        v
    end
end

function display_results(df)
    case_starts = [findfirst(isequal(c), collect(zip(df.case, df.repeat)))
                   for c in unique(collect(zip(df.case, df.repeat)))]

    function is_maxcase(df, columns)
        function inner(data, i, j)
            !(j in columns) && return false
            ismissing(data[i, j]) && return false
            maxval = maximum(filter(d -> !ismissing(d), df[(df.case .== data[i, 1]) .& (df.repeat .== data[i, 2]), j]))
            data[i, j] == maxval
        end
    end
    function is_mincase(df, columns)
        function inner(data, i, j)
            !(j in columns) && return false
            ismissing(data[i, j]) && return false
            minval = minimum(filter(d -> !ismissing(d), df[(df.case .== data[i, 1]) .& (df.repeat .== data[i, 2]), j]))
            data[i, j] == minval
        end
    end

    hlmaxcond = Highlighter(f=is_maxcase(df, [6, 7]), crayon=crayon"red")
    hlmincond = Highlighter(f=is_mincase(df, [6, 7]), crayon=crayon"blue")
    pretty_table(df,
                 body_hlines=case_starts[2:end] .- 1,
                 formatters=(ft_sigdigits([3, 7, 8, 9]), ft_iterations([6])),
                 highlighters=(hlmaxcond, hlmincond),
                 nosubheader=true,
                )
end


function prepare_table(df; cases, mixings, latex=false, shorten_resta=false)
    header = ["case", "repeat"]
    subheader = ["", ""]
    for mix in mixings
        shorten_resta && mix == "RestaGaAs" && (mix = "Resta")
        mixing_clean = get(mixing_cleanup, mix, mix)
        append!(header, [mixing_clean, mixing_clean])
        append!(subheader, ["it", "κ"])
    end
    header = permutedims(hcat(header, subheader), (2, 1))

    rows = []
    itcols = Int[]
    condcols = Int[]
    caserows = Int[]
    df_filter = df[in.(df.case, Ref(cases)), :]
    for case in cases
        df_case = df_filter[df_filter.case .== case, :]
        first = true
        for repeat in unique(df_case.repeat)
            df_repeat = df_case[df_case.repeat .== repeat, :]

            row = [first ? get(case_cleanup, case, case) : "", repeat]
            first = false
            for mix in mixings
                df_mix = df_repeat[df_repeat.mixing .== mix, :]
                push!(itcols,   1 + length(row))
                push!(condcols, 2 + length(row))
                if isempty(df_mix)
                    append!(row, Any[missing, missing])
                else
                    @assert size(df_mix, 1) == 1
                    append!(row, Any[df_mix.iterations[1], df_mix.condition_number[1]])
                end
            end
            push!(rows, row)
        end
        push!(caserows, length(rows))
    end
    itcols = unique(itcols)
    condcols = unique(condcols)

    data = permutedims(hcat(rows...), (2, 1))
    kwargs = (
        formatters=(ft_sigdigits(condcols), ft_iterations(itcols)),
        vlines=itcols .- 1,
    )
    if !latex
        hlcond = Highlighter((data,i,j)->(j in condcols) && !ismissing(data[i, j]),
                             (h,data,i,j)->begin
                                 condrange = (minimum([d for d in data[i, condcols]
                                                       if !ismissing(d)]),
                                              maximum([d for d in data[i, condcols]
                                                       if !ismissing(d)]))
                                 color = get(ColorSchemes.colorschemes[:coolwarm],
                                              data[i,j], condrange)
                                 Crayon(foreground=(round(Int, color.r * 255),
                                                    round(Int, color.g * 255),
                                                    round(Int, color.b * 255)))
                              end)
        hlit = Highlighter((data,i,j)->(j in itcols) && !ismissing(data[i, j]),
                           (h,data,i,j)->begin
                               itrange = (minimum([d for d in data[i, itcols]
                                                   if !ismissing(d)]), 50)
                               color = get(ColorSchemes.colorschemes[:coolwarm],
                                           data[i,j], itrange)
                               Crayon(foreground=(round(Int, color.r * 255),
                                                  round(Int, color.g * 255),
                                                  round(Int, color.b * 255)))
                            end)
        pretty_table(data, header; highlighters=(hlcond, hlit),
                     hlines=append!([1], caserows[1:end-1] .+ 1),
                     kwargs...)
    else
        function is_maxcase(columns)
            function inner(data, i, j)
                !(j in columns) && return false
                ismissing(data[i, j]) && return false
                data[i, j] == maximum([d for d in data[i, columns] if !ismissing(d)])
            end
        end
        function is_mincase(columns)
            function inner(data, i, j)
                !(j in columns) && return false
                ismissing(data[i, j]) && return false
                data[i, j] == minimum([d for d in data[i, columns] if !ismissing(d)])
            end
        end

        tf = LatexTableFormat(top_line="",
                              #header_line=raw"\midrule",
                              #mid_line="\n" * raw"\midrule",
                              bottom_line="",
                              header_envs=[],
                              subheader_envs=[],
                             )
        highlighters = (
            LatexHighlighter(is_maxcase(condcols), ["color{red}"]),
            LatexHighlighter(is_mincase(condcols), ["color{blue}"]),
            LatexHighlighter(is_maxcase(itcols), ["color{red}"]),
            LatexHighlighter(is_mincase(itcols), ["color{blue}"]),
        )
        out = pretty_table(String, data, header; backend=:latex, highlighters=highlighters,
                           hlines=caserows[1:end-1],
                           tf=tf, kwargs...)
        for (i, mix) in enumerate(mixings)
            shorten_resta && mix == "RestaGaAs" && (mix = "Resta")
            mc = get(mixing_cleanup, mix, mix)
            bar = i == length(mixings) ? "" : "|"
            out = replace(out, "$mc & $mc" => "\\multicolumn{2}{c$bar}{$mc}")
        end
        out = replace(out, "κ"  => raw"$\kappa$")
        out = replace(out, "SiO2"  => raw"\ce{SiO2}")
        out = replace(out, ">50"   => raw"$>50$")
        out = replace(out, "missing"  => raw"\textcolor{green}{\textbf{??}}")

        out
    end
end


function table_bulk(df; kwargs...)
    prepare_table(df; cases=["SiO2", "GaAs", "Al"],
                  mixings=["None", "RestaSiO2", "RestaGaAs", "Kerker"], kwargs...)
end


function table_mixed(df; kwargs...)
    prepare_table(df; cases=["SiO2HVac", "GaAsVac", "AlVac", "GaAsSiO2H", "AlSiO2H",
                             "AlGaAs", "AlGaAsSiO2H"],
                  mixings=["None", "RestaGaAs", "Kerker", "HybridVac",
                           "HybridGaAs", "Custom"], shorten_resta=true, kwargs...)
end


function plot_scatter_conditioning(df)
    bulk = ["SiO2", "GaAs", "Al"]
    mixed = ["SiO2HVac", "GaAsVac", "AlVac", "GaAsSiO2H", "AlSiO2H",
             "AlGaAs", "AlGaAsSiO2H"]

    df = df[(.!ismissing.(df.iterations)) .& (df.iterations .!= typemax(Int)), :]
    df_bulk = df[in.(df.case, Ref(bulk)), :]
    df_mixed = df[in.(df.case, Ref(mixed)), :]

    p = plot(legend=:topleft)
    scatter!(p, df_bulk.condition_number, df_bulk.iterations, m=:x, label="bulk")
    scatter!(p, df_mixed.condition_number, df_mixed.iterations, m=:x, label="mixed")

    # Theoretical convergence rate for arnoldi methods
    conv_rate(κ) = (sqrt(κ) - 1) / (sqrt(κ) + 1)
    n_iter(κ) = -11 / log(conv_rate(κ))

    κs = collect(1:0.1:160)
    plot!(κs, n_iter.(κs), label="-11 / log(rate(κ))", ls=:dot)

    p
end
