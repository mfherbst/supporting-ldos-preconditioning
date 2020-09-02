config[@__DIR__] = (
    directory=@__DIR__,
    case="GaAs",
    n_repeat=20,
    kwargs_scf=(randomise_guess=true, ),
    firstmixkey="RestaGaAs",
    #
    construct_localiser=(lattice, atoms) -> identity,
    mixings=Dict(
        "Kerker"     => (α, localiser) -> KerkerMixing(α=α, kTF=1.0),
        "None"       => (α, localiser) -> SimpleMixing(α=α),
        "RestaGaAs"  => (α, localiser) -> RestaMixing( α=α, kTF=1.0, εr=14),
        "RestaSiO2"  => (α, localiser) -> RestaMixing( α=α, kTF=1.0, εr=2.4),
    ),
)
