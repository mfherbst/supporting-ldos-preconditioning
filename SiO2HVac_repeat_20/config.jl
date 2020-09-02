import SpecialFunctions: erf

function construct_localiser_SiO2HVac(lattice, atoms; σ=30)
    @assert lattice == Diagonal(diag(lattice))
    @assert maximum(lattice) == lattice[3, 3]  # 3rd axis is the elongated one
    @assert length(atoms) == 3                 # Two types of atoms
    @assert atoms[1][1].symbol == :Si          # Silicon symbol
    @assert atoms[2][1].symbol == :O           # Oxygen symbol
    @assert atoms[3][1].symbol == :H           # Hydrogen symbol
    @assert σ ≥ 15

    # The slabs are usually from 0.0 through to around 0.5
    # Using this we find first and last silicon atom
    last_z  = sort([r[3] for iat in (1, 2, 3) for r in atoms[iat][2]], by=z -> abs(z - 0.5))[1]
    first_z = sort([r[3] for iat in (1, 2, 3) for r in atoms[iat][2]],
                   by=z -> min(abs(z - 0.0), abs(z - 1.0)))[1]
    first_z > 0.75 && (first_z -= 1.0)

    # https://en.wikipedia.org/wiki/Atomic_radii_of_the_elements_%28data_page%29
    atomic_radii = Dict(:O => 0.48 * DFTK.units.Å, :Si => 1.11 * DFTK.units.Å, :H => 0.53)
    max_at = max(atomic_radii[:O], atomic_radii[:Si], atomic_radii[:H])
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
    case="SiO2HVac",
    n_repeat=20,
    kwargs_scf=(kspacing=0.15, tol=5e-11, ),
    firstmixkey="Custom",
    #
    construct_localiser=construct_localiser_SiO2HVac,
    mixings=Dict(
        "Custom"      => (α, localiser) -> HybridMixing(α=α, kTF=1.0, εr=14, localisation=localiser),
        "HybridGaAs"  => (α, localiser) -> HybridMixing(α=α, kTF=1.0, εr=14, localisation=identity),
        "HybridVac"   => (α, localiser) -> HybridMixing(α=α, kTF=1.0, εr=1,  localisation=identity),
        "Kerker"      => (α, localiser) -> KerkerMixing(α=α, kTF=1.0),
        "None"        => (α, localiser) -> SimpleMixing(α=α),
        "RestaGaAs"   => (α, localiser) -> RestaMixing( α=α, kTF=1.0, εr=14),
    ),
)
