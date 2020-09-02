config[@__DIR__] = (
    directory=@__DIR__,
    case="AlVac",
    n_repeat=10,
    kwargs_scf=(kspacing=0.15, smearing=Smearing.Gaussian(), temperature=0.001,
                determine_bands=nelec -> nelec // 2 + 2sqrt(nelec // 2), ),
    firstmixkey="HybridVac",
    #
    construct_localiser=(lattice, atoms) -> identity,
    mixings=Dict(
        # "Custom"      => same as HybridVac
        "HybridGaAs"  => (α, localiser) -> HybridMixing(α=α, kTF=1.0, εr=14, localisation=identity),
        "HybridVac"   => (α, localiser) -> HybridMixing(α=α, kTF=1.0, εr=1,  localisation=identity),
        "Kerker"      => (α, localiser) -> KerkerMixing(α=α, kTF=1.0),
        "None"        => (α, localiser) -> SimpleMixing(α=α),
        "RestaGaAs"   => (α, localiser) -> RestaMixing( α=α, kTF=1.0, εr=14),
    ),
)
