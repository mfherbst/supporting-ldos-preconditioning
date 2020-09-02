config[@__DIR__] = (
    directory=@__DIR__,
    case="GaAs",
    n_repeat=40,
    kwargs_scf=(randomise_guess=true, random_ampl=1e-3, ),
    firstmixkey="RestaGaAs",
    # firstα=0.1,
    #
    construct_localiser=(lattice, atoms) -> identity,
    mixings=Dict(
        "Kerker"     => (α, localiser) -> KerkerMixing(α=α, kTF=1.0),
        "None"       => (α, localiser) -> SimpleMixing(α=α),
        "RestaGaAs"  => (α, localiser) -> RestaMixing( α=α, kTF=1.0, εr=14),
        "RestaSiO2"  => (α, localiser) -> RestaMixing( α=α, kTF=1.0, εr=2.4),
    ),
)
