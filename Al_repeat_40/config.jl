config[@__DIR__] = (
    directory=@__DIR__,
    case="Al",
    n_repeat=40,
    kwargs_scf=(smearing=Smearing.Gaussian(), temperature=0.001, randomise_guess=true),
    firstmixkey="Kerker",
    #
    construct_localiser=(lattice, atoms) -> identity,
    mixings=Dict(
        "Kerker"     => (α, localiser) -> KerkerMixing(α=α, kTF=1.0),
        "None"       => (α, localiser) -> SimpleMixing(α=α),
        "RestaGaAs"  => (α, localiser) -> RestaMixing( α=α, kTF=1.0, εr=14),
        "RestaSiO2"  => (α, localiser) -> RestaMixing( α=α, kTF=1.0, εr=2.4),
    ),
)
