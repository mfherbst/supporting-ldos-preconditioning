import SpecialFunctions: erf

# Construct custom localiser for the GaAs area
function construct_localiser_GaAsSiO2H(lattice, atoms; σ=30)
    @assert lattice == Diagonal(diag(lattice))
    @assert maximum(lattice) == lattice[3, 3]  # 3rd axis is the elongated one
    @assert length(atoms) == 5                 # Two types of atoms
    @assert atoms[1][1].symbol == :Ga          # Gallium symbol
    @assert atoms[2][1].symbol == :As          # Arsenic symbol
    @assert σ ≥ 15

    # The slabs are usually from 0.0 through to around 0.5
    # Using this we find first and last silicon atom
    last_z  = sort([r[3] for iat in (1, 2) for r in atoms[iat][2]], by=z -> abs(z - 0.5))[1]
    first_z = sort([r[3] for iat in (1, 2) for r in atoms[iat][2]],
                   by=z -> min(abs(z - 0.0), abs(z - 1.0)))[1]
    first_z > 0.75 && (first_z -= 1.0)


    # https://en.wikipedia.org/wiki/Atomic_radii_of_the_elements_%28data_page%29
    atomic_radii = Dict(:Ga => 1.36 * DFTK.units.Å, :As => 1.14 * DFTK.units.Å)
    max_at = max(atomic_radii[:Ga], atomic_radii[:As])
    last_z  += max_at / lattice[3, 3]  # Convert to fractional coords
    first_z -= max_at / lattice[3, 3]

    @assert -0.25 < first_z < 0.25
    @assert 0.25 < last_z < 0.75
    periodise(x) = ifelse(x < 0.75, x + 1, x)

    first_z = periodise(first_z)
    last_z = periodise(last_z)
    r -> max(0, + 0.5 * erf(σ * (periodise(r[3]) - first_z))
                - 0.5 * erf(σ * (periodise(r[3]) - last_z)))
end


config[@__DIR__] = (
    directory=@__DIR__,
    case="GaAsSiO2H",
    n_repeat=20,
    kwargs_scf=(kspacing=0.15, tol=5e-11, ),
    firstmixkey="Custom",
    #
    construct_localiser=construct_localiser_GaAsSiO2H,
    mixings=Dict(
        "Custom"      => (α, localiser) -> HybridMixing(α=α, kTF=1.0, εr=14, localisation=localiser),
        "HybridGaAs"  => (α, localiser) -> HybridMixing(α=α, kTF=1.0, εr=14, localisation=identity),
        "HybridVac"   => (α, localiser) -> HybridMixing(α=α, kTF=1.0, εr=1,  localisation=identity),
        "Kerker"      => (α, localiser) -> KerkerMixing(α=α, kTF=1.0),
        "None"        => (α, localiser) -> SimpleMixing(α=α),
        "RestaGaAs"   => (α, localiser) -> RestaMixing( α=α, kTF=1.0, εr=14),
    ),
)
