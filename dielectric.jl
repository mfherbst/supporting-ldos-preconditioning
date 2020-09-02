# Compute a few eigenvalues of the dielectric matrix at (q=0, ω=0) iteratively
using DFTK
using Plots
import KrylovKit: ArnoldiIterator, Orthogonalizer, OrthonormalBasis, KrylovDefaults, orthogonalize!
using KrylovKit
using Dates
using HDF5
using WriteVTK
using Printf

#
# For next time: Consider residual error in stopping criterion and
# stop only if estimated error in κ is less than 5%
#

# KrylovKit custom orthogonaliser
struct OrthogonalizeAndProject{F, O <: Orthogonalizer} <: Orthogonalizer
    projector!::F
    orth::O
end
OrthogonalizeAndProject(projector) = OrthogonalizeAndProject(projector, KrylovDefaults.orth)
function KrylovKit.orthogonalize!(v::T, b::OrthonormalBasis{T}, x::AbstractVector, alg::OrthogonalizeAndProject) where {T}
    v, x = orthogonalize!(v, b, x, alg.orth)
    v .= alg.projector!(v)::T
    v, x
end
function KrylovKit.orthogonalize!(v::T, q::T, alg::OrthogonalizeAndProject) where {T}
    v, s = orthogonalize!(v, q, alg.orth)
    v .= alg.projector!(v)::T
    v, s
end


# A straightfoward Arnoldi eigensolver that diagonalizes the matrix at each step
# This is more efficient than Arpack when `f` is very expensive
# Convergence is flagged after a given number of iterations or if the condition number
# changes by less than a given relative tolerance
function arnoldi(f, x0; cond_rtol=1e-3, iterations=30, min_iterations=10, n_print=5,
                 projector=identity)
    println()
    println("Starting Arnoldi (now = $(Dates.now()) ...")
    println()
    flush(stdout)

    oldcond = NaN
    for (V, B, r, nr, b) in ArnoldiIterator(f, x0, OrthogonalizeAndProject(projector))
        @assert isreal(r)

        # A * V = V * B + r * b'
        V = hcat(V...)
        AV = V*B + r*b'

        ew, ev = eigen(B, sortby=real)
        Vr = V*ev
        AVr = AV*ev
        R = AVr - Vr * Diagonal(ew)

        N = size(V, 2)
        normr = [norm(r) for r in eachcol(R)]

        println("#--- $N ---#")
        println("idcs      evals   residnorms   |P(X)|")
        inds = unique(append!(collect(1:min(N, n_print)), max(1, N-n_print):N))
        for i in inds
            Xi = @view Vr[:, i]
            norm_asym = abs(imag(ew[i])) < 1e-6 ? norm(Xi - projector(Xi)) : NaN
            @printf "% 3i  %10.6g  %10.6g  %10.6g\n" i real(ew[i]) normr[i] norm_asym
        end
        any(imag.(ew[inds]) .> 1e-5) && println("Warn: Suppressed imaginary part.")
        println()
        flush(stdout)

        cond = real(ew[end] / ew[1])
        is_converged = abs(cond - oldcond) / oldcond < cond_rtol && N > min_iterations
        oldcond = cond
        if N ≥ iterations || is_converged
            return (λ=ew, X=Vr, AX=AVr, residual_norms=normr)
        end
    end
end


function dump_arnoldi_hdf5(scfres, arnoldires, mixing, filename)
    gridsize = size(scfres.ρ.real)
    n_epairs = length(arnoldires.λ)
    isfile(filename) && rm(filename)
    h5open(filename, "w") do h5f
        arnoldi = g_create(h5f, "arnoldi")
        arnoldi["λ"]               = arnoldires.λ
        arnoldi["residual_norms"]  = arnoldires.residual_norms
        arnoldi["X",  "shuffle", (), "compress", 3] = real(reshape(arnoldires.X,
                                                                   gridsize..., n_epairs))
        arnoldi["AX", "shuffle", (), "compress", 3] = real(reshape(arnoldires.AX,
                                                                   gridsize..., n_epairs))
    end
end


function dump_arnoldi_vtk(scfres, arnoldires, mixing, filename)
    gridsize = size(scfres.ρ.real)
    arrays = Array{Float64, 3}[]
    labels = String[]

    for i in 1:length(arnoldires.λ)
        push!(arrays, real(reshape(arnoldires.X[:, i], gridsize)))
        push!(labels, "X_$(@sprintf "%02d" i)")

        push!(arrays, real(reshape(arnoldires.AX[:, i], gridsize)))
        push!(labels, "AX_$(@sprintf "%02d" i)")
    end

    # Save the LDOS
    if scfres.ham.basis.model.temperature > 0
        ldos = LDOS(scfres.εF, scfres.ham.basis, scfres.eigenvalues, scfres.ψ)
        push!(arrays, ldos)
        push!(labels, "LDOS")
    end

    idx = findfirst(isa.(scfres.ham.basis.model.term_types, AtomicLocal))
    push!(arrays, scfres.ham.basis.terms[idx].potential)
    push!(labels, "pot_atomic_local")

    push!(arrays, DFTK.total_local_potential(scfres.ham))
    push!(labels, "pot_total_local")

    push!(arrays, scfres.ρ.real)
    push!(labels, "ρ")
    vtk_write_array(filename, Tuple(arrays), Tuple(labels))
end


function run_arnoldi_dielectric(scfres; mixing=nothing, iterations=30, cond_rtol=5e-4,
                                fileprefix="arnoldi", force_lda=false)
    if isfile(fileprefix * ".hdf5")
        arnoldifile = fileprefix * ".hdf5"
        return (λ=h5read(arnoldifile, "arnoldi/λ"),
                residual_norms=h5read(arnoldifile, "arnoldi/residual_norms"))
    end

    # self_consistent_field uses 3 extra bands, which are not converged by the eigensolver
    # For the χ0 application, bands need to be perfectly converged, so we discard them here
    ψ_cvg = [@view ψk[:, 1:end-3]  for ψk in scfres.ψ]
    eigenvalues_cvg = [εk[1:end-3] for εk in scfres.eigenvalues]

    ham = scfres.ham
    basis = ham.basis
    if force_lda
        origmodel = basis.model
        model = model_LDA(origmodel.lattice, origmodel.atoms, smearing=origmodel.smearing,
                          temperature=origmodel.temperature)
        kcoords = [kpt.coordinate for kpt in basis.kpoints]
        basis = PlaneWaveBasis(model, basis.Ecut, kcoords, basis.ksymops, basis.symops,
                               fft_size=basis.fft_size)
    end

    # Apply ε = 1 - χ0 (vc + fxc)
    function eps_fun(dρ)
        dρ = reshape(dρ, size(scfres.ρ.real))
        dρ = from_real(basis, dρ)
        dv = apply_kernel(basis, dρ; ρ=scfres.ρ)
        χdv = apply_χ0(ham, ψ_cvg, scfres.εF, eigenvalues_cvg, dv)
        vec((dρ - χdv).real)
    end
    function pinv_fun(dρ)
        dρ = reshape(dρ, size(scfres.ρ.real))
        dρ_dc = sum(dρ) / length(dρ)
        dρ .-= dρ_dc

        # Make dummy ρin and ρout to use mixing interface,
        # such that δρ = ρout - ρin is as required
        ρin  = scfres.ρ
        ρout = from_real(basis, dρ .+ scfres.ρ.real)
        ρnext = DFTK.mix(mixing, basis, ρin, ρout; scfres...)

        ρdiff = ρnext - ρin
        vec(ρdiff.real .+ dρ_dc)
    end

    # Projector into the subspace of interest
    # (i.e. modes symmetric wrt. the lattice symmetries of the problem)
    function projector(x)
        x_dc = sum(x) / length(x)
        x .-= x_dc
        x_real = reshape(real(x), size(scfres.ρ.real))
        vec(DFTK.symmetrize(from_real(basis, x_real)).real)
    end

    apply = projector ∘ eps_fun
    if !isnothing(mixing)
        apply = projector ∘ pinv_fun ∘ eps_fun
    end

    # Construct guess: Normalised, mostly random, but large mass on smallest G mode
    aguess = 0.1randn(size(scfres.ρ.real))
    Gmin = minimum(G -> iszero(G) ? Inf : norm(model.recip_lattice * G), G_vectors(basis))
    idx_Gmin = findfirst(G -> abs(norm(model.recip_lattice * G) - Gmin) < 1e-12,
                         collect(G_vectors(basis)))
    aguess[idx_Gmin] = 2.0
    aguess /= norm(aguess)
    aguess = projector(aguess)
    aguess = projector(aguess)

    arnoldires = arnoldi(apply, vec(aguess), cond_rtol=cond_rtol, iterations=iterations,
                         projector=projector)
    dump_arnoldi_hdf5(scfres, arnoldires, mixing, fileprefix * ".hdf5")
    dump_arnoldi_vtk(scfres, arnoldires, mixing, fileprefix)

    arnoldires
end
