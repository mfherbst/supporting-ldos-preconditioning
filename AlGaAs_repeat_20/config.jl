config[@__DIR__] = (
    directory=@__DIR__,
    case="AlGaAs",
    n_repeat=20,
    kwargs_scf=(kspacing=0.15, smearing=Smearing.Gaussian(), temperature=0.001),
    firstmixkey="HybridGaAs",
    #
    construct_localiser=(lattice, atoms) -> identity,
    mixings=Dict(
        # "Custom"      => same as HybridGaAs
        "HybridGaAs"  => (α, localiser) -> HybridMixing(α=α, kTF=1.0, εr=14, localisation=identity),
        "HybridVac"   => (α, localiser) -> HybridMixing(α=α, kTF=1.0, εr=1,  localisation=identity),
        "Kerker"      => (α, localiser) -> KerkerMixing(α=α, kTF=1.0),
        "None"        => (α, localiser) -> SimpleMixing(α=α),
        "RestaGaAs"   => (α, localiser) -> RestaMixing( α=α, kTF=1.0, εr=14),
    ),
)
