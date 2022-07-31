using ExponentialUtilities: arnoldi!, KrylovSubspace, expv!, phiv!
using FastExpm: fastExpm
using ReachabilityAnalysis: reach_homog_LGG09!, step_size, zeroI, ReachSolution, _Eplus

# Computes: Ω0 = CH(X0, ΦX0 ⊕ Eplus(X0, δ))
# using Krylov for Eplus and FastExpm for Phi
# Currently it is assuming that X0 is a singleton,
# but we should generalize it
function solve_krylov_discr(prob; NSTEPS, alg::LGG09{N, AM, VN, TN}) where {N, AM, VN, TN}
    X0 = if prob.x0 isa AbstractHyperrectangle
        prob.x0
    else
        box_approximation(prob.x0)
    end
    s = prob.s
    # Currently only "C" case is supported.
    @assert s isa LinearContinuousSystem
    A = sparse(s.A)
    δ = step_size(alg)

    # discretization
    Ep = _Eplus(A, X0, δ)
    Φ = real.(fastExpm(A .* δ, threshold=eps(Float64), nonzero_tol=eps(Float64)))
    Ω0 = CH(X0, Φ*X0 ⊕ Ep)

    # preallocate output flowpipe
    SN = SubArray{N, 1, Matrix{N}, Tuple{Base.Slice{Base.OneTo{Int}}, Int}, true}
    RT = TemplateReachSet{N, VN, TN, SN}
    F = Vector{RT}(undef, NSTEPS)

    # propagate
    X = Universe(dim(X0))
    ρℓ = reach_homog_LGG09!(F, alg.template, Ω0, Φ, NSTEPS, δ, X, zeroI, Val(alg.cache), Val(alg.threaded))

    # output solution
    dict = Dict{Symbol, Any}(:sfmat => ρℓ, :alg_vars => alg.vars)
    ReachSolution(Flowpipe(F, dict), alg)
end
