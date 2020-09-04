# Black-box inhomogeneous preconditioning for self-consistent field iterations in density-functional theory

[![](https://img.shields.io/badge/arxiv-2009.01665-red)](https://arxiv.org/abs/2009.01665)

Supporting information containing structures,
raw data and computational scripts for the paper:

Michael F. Herbst and Antoine Levitt  
*Black-box inhomogeneous preconditioning for self-consistent field iterations in density-functional theory*  
Preprint on [arxiv (2009.01665)](https://arxiv.org/abs/2009.01665)

The code in this repository has been used to run all calculations
and produce all plots of the above paper.
It relies on [DFTK](https://dftk.org) version 0.1.8
and [ASE](https://wiki.fysik.dtu.dk/ase/).
For more details see
[this blog article](https://michael-herbst.com/2020-ldos-preconditioning.html).

In particular:
  - The folder `systems` contains the employed lattice setups.
    See the `load_system` function of [builder.jl](builder.jl) for a parser to DFTK
    data structures.
  - `<case>/config.jl` contain the parameters for each test case, [dftk.jl](dftk.jl)
    the general DFTK parameters.
  - [dielectric.jl](dielectric.jl) contains the details of our Arnoldi procedure
    to compute the eigenvalues of the dielectric operator.
  - [plot.jl](plot.jl) contains the translation from the "internal" preconditioner names
    (used for the scripts in the repository) and the "published" names (used in the paper).

## Running the code and reproducing the plots
Running the code requires an installation of
[Julia 1.4.0](https://julialang.org/downloads/#current_stable_release),
of [DFTK](https://docs.dftk.org/dev/guide/installation/)
and of [ASE](https://wiki.fysik.dtu.dk/ase/).
With this setup can generate the plots of the paper by executing:
```bash
julia --project=@. -e "import Pkg; Pkg.instantiate()"  # Install dependencies
julia run.jl   # Generate data
julia plot.jl  # Generate plots
```

Be aware that generating the data takes a long time
(like two weeks on a couple of cluster nodes). Most most raw data is,
however, included in the repository, such that some of the plotting works
without running the calculations beforehand.

### Quantum ESPRESSO setup
For a small number of cases (see the paper) we compare against
[Quantum ESPRESSO](http://www.quantum-espresso.org/).
Running the respective `<case>/run_qe.jl` scripts requires
a working Quantum ESPRESSO setup.
We use the program via the ASE calculator interface.
For this interface to find Quantum ESPRESSO you will need to adapt
[quantum_espresso.jl](quantum_espresso.jl) to suit your Quantum ESPRESSO
installation.
