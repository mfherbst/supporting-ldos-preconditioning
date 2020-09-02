include("../dftk.jl")
include("../builder.jl")
import FFTW
using LinearAlgebra
import TimerOutputs

function run(stepsize::Integer)
    FFTW.set_num_threads(4)
    BLAS.set_num_threads(4)

    lattice, atoms = load_system("AlVac", 10)
    mixing         = KerkerMixing(α=1.84, kTF=1.0)  # optimal damping for Kerker
    prefix         = joinpath(@__DIR__, "anderson_$stepsize")

    logfile = joinpath(@__DIR__, "anderson_$stepsize.log")
    println("Anderson m=$stepsize  ->  $logfile")
    open(logfile, "w") do fp
        redirect_stdout(fp) do
            println("#\n#-- Anderson m=$stepsize\n#")
            TimerOutputs.reset_timer!(DFTK.timer)
            run_dftk(lattice, atoms; fileprefix=prefix, mixing=mixing,
                     smearing=Smearing.Gaussian(), temperature=0.001,
                     kspacing=0.15, 
                     randomise_guess=false, stepsize=stepsize, maxiter=100)
            println(DFTK.timer)
        end
    end
end
