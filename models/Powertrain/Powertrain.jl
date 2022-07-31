using ReachabilityAnalysis, SparseArrays
using ReachabilityAnalysis: add_dimension

LazySets.set_ztol(Float64, 1e-14)

using LazySets: HalfSpace

# physical units of the state variables
# [x₁] = rad
# [x₂] = Nm
# [x₃] = rad
# [x₄] = rad/s
# [x₅] = rad
# [x₆] = rad/s
# [x₇] = rad/s
# [x₈] = rad
# [x₉] = rad/s
# ...
# [x_(2θ+6)] = rad
# [x_(2θ+7)] = rad/s

# physical constants
# indices 'm' (ₘ) and 'l' (ₗ) refer to the motor and the load
# indices 'i' (ᵢ) refer to the numbering of additional rotating masses

const t_init = 0.2  # time to stay in the initial location
const α = 0.03  # backlash size (half of the gap width)
const τ_eng = 0.1  # engine time constant
const γ = 12.0 # gearbox ratio (dimensionless)
const u = 5.0 # requested engine torque

# moments of inertia [kg m²]
const Jₗ = 140.
const Jₘ = 0.3
const Jᵢ = 100.0  # TODO the paper says 0.01, the SpaceEx model uses 100

# viscous friction constants [Nm s/rad]
const bₗ = 5.6
const bₘ = 0.0
const bᵢ = 1.0

# shaft stiffness [Nm/rad]
const kᵢ = 1e5
const kₛ = 1e4

# PID parameters
const k_P = 0.5  # [Nms/rad]
const k_I = 0.5  # [Nm/rad]
const k_D = 0.5  # [Nms²/rad]

function print_dynamics(A, b, location_name)
    println("dynamics of location $location_name:")
    for i in 1:size(A, 1)-1  # ignore the last dimension (time)
        print("x_$i' = ")
        for j in 1:size(A, 2)
            if !iszero(A[i,j])
                print("$(A[i,j]) x_$j + ")
            end
        end
        println("$(b[i])\n")
    end
end

function get_dynamics(kₛ, α, u, n, θ)
    # linear dynamics
    A = spzeros(n, n)

    A[1, 7] = 1.0 / γ
    A[1, 9] = -1.0

    A[2, 1] = (-k_I * γ + k_D * kₛ / (γ * Jₘ)) / τ_eng
    A[2, 2] = (-k_D / Jₘ - 1.0) / τ_eng
    A[2, 3] = k_I * γ / τ_eng
    A[2, 4] = k_P * γ / τ_eng
    A[2, 7] = (-k_P + k_D * bₘ /Jₘ) / τ_eng
    A[2, 8] = -k_I * γ / τ_eng

    A[3, 4] = 1.0

    A[5, 6] = 1.0

    A[6, 5] = -kᵢ / Jₗ
    A[6, 6] = -bₗ / Jₗ
    A[6, 2*θ+6] = kᵢ / Jₗ

    A[7, 1] = -kₛ / (Jₘ * γ)
    A[7, 2] = 1.0 / Jₘ
    A[7, 7] = -bₘ / Jₘ

    i = 8
    while i < n-1
        A[i, i+1] = 1.0

        if i == 8
            # x9 has special dynamics
            A[i+1, 1] = kₛ / Jᵢ
            A[i+1, i] = -kᵢ / Jᵢ
        else
            A[i+1, i-2] = kᵢ / Jᵢ
            A[i+1, i] = -2. * kᵢ / Jᵢ
        end
        A[i+1, i+1] = -bᵢ / Jᵢ

        # wrap-around to x5 in the last step
        j = (i == n-2) ? 5 : i+2
        A[i+1, j] = kᵢ / Jᵢ

        i += 2
    end

    # affine vector
    b = spzeros(n)
    b[2] = k_D * (γ * u - kₛ * α / (Jₘ * γ)) / τ_eng
    b[4] = u
    b[7] = kₛ * α / (Jₘ * γ)
    b[9] = -kₛ * α / Jᵢ
    b[n] = 1.0  # time

    return A, b
end

function get_initial_condition(n, X0_scale)
    c = Vector{Float64}(undef, n)
    g = Vector{Float64}(undef, n)
    c[1:7] = [-0.0432, -11., 0., 30., 0., 30., 360.]
    g[1:7] = [0.0056, 4.67, 0., 10., 0., 10., 120.]
    i = 8
    while i < n
        c[i] = -0.0013
        g[i] = 0.0006
        i += 1
        c[i] = 30.
        g[i] = 10.
        i += 1
    end
    c[n] = 0.0
    g[n] = 0.0
    if X0_scale < 1.0
        g = X0_scale * g
    end

    X0 = Zonotope(c, hcat(g))
    return X0
end

function powertrain(; θ::Int=1, X0_scale::Float64=1.0)
    @assert θ > 0 "θ must be positive, but was $θ"
    @assert (X0_scale > 0.0 && X0_scale <= 1.0) "scale $X0_scale ∉ (0, 1]"

    # dimension of state space (last dimension is time)
    n = 2 * θ + 7 + 1

    # hybrid automaton
    automaton = GraphAutomaton(4)

    # negAngle
    A, b = get_dynamics(kₛ, -α, u, n, θ)
    X = HalfSpace(sparsevec([1], [1.], n), -α)  # x1 <= -α
    m_negAngle = @system(x' = A * x + b, x ∈ X)

    # deadzone
    A, b = get_dynamics(0., -α, u, n, θ)
    X = HPolyhedron([HalfSpace(sparsevec([1], [-1.], n), α),  # x1 >= -α
                     HalfSpace(sparsevec([1], [1.], n), α)])  # x1 <= α
    m_deadzone = @system(x' = A * x + b, x ∈ X)

    # posAngle
    A, b = get_dynamics(kₛ, α, u, n, θ)
    X = HalfSpace(sparsevec([1], [-1.], n), -α)  # x1 >= α
    m_posAngle = @system(x' = A * x + b, x ∈ X)

    # negAngleInit
    A, b = get_dynamics(kₛ, -α, -u, n, θ)
    X = HalfSpace(sparsevec([n], [1.], n), t_init)  # t <= t_init
    m_negAngleInit = @system(x' = A * x + b, x ∈ X)

    # modes
    modes = [m_negAngle, m_deadzone, m_posAngle, m_negAngleInit]

    # transition negAngleInit -> negAngle
    add_transition!(automaton, 4, 1, 1)
    guard = HalfSpace(sparsevec([n], [-1.], n), -t_init)  # t >= t_init
    r_41 = ConstrainedIdentityMap(n, guard)

    # transition negAngle -> deadzone
    add_transition!(automaton, 1, 2, 2)
    guard = HalfSpace(sparsevec([1], [-1.], n), α)  # x1 >= -α
    r_12 = ConstrainedIdentityMap(n, guard)

    # transition deadzone -> posAngle
    add_transition!(automaton, 2, 3, 3)
    guard = HalfSpace(sparsevec([1], [-1.], n), -α)  # x1 >= α
    r_23 = ConstrainedIdentityMap(n, guard)

    # TODO the SpaceEx model does not contain the following transitions
#     # transition deadzone -> negAngle
#     add_transition!(automaton, 2, 1, 4)
#     guard = HalfSpace(sparsevec([1], [1.], n), -α)  # x1 <= -α
#     r_21 = ConstrainedIdentityMap(n, guard)
#     # transition posAngle -> deadzone
#     add_transition!(automaton, 3, 2, 5)
#     guard = HalfSpace(sparsevec([1], [1.], n), α)  # x1 <= α
#     r_32 = ConstrainedIdentityMap(n, guard)

    # transition annotations
    resetmaps = [r_41, r_12, r_23]

    # switching
    switchings = [HybridSystems.AutonomousSwitching()]

    H = HybridSystem(automaton, modes, resetmaps, switchings)

    # initial condition TODO check initial mode
    X0 = get_initial_condition(n, X0_scale)
    initial_condition = [(4, X0)]

    return InitialValueProblem(H, initial_condition)
end

function powertrain_homog(; θ::Int=1, X0_scale::Float64=1.0)
    @assert θ > 0 "θ must be positive, but was $θ"
    @assert (X0_scale > 0.0 && X0_scale <= 1.0) "scale $X0_scale ∉ (0, 1]"

    # dimension of state space (last dimension is time)
    n = 2 * θ + 7 + 1

    # hybrid automaton
    automaton = GraphAutomaton(4)

    # negAngle
    A, b = get_dynamics(kₛ, -α, u, n, θ)
    X = HalfSpace(sparsevec([1], [1.], n+1), -α)  # x1 <= -α
    Aext = add_dimension(A, 1)
    Aext[1:n, n+1] .= b
    m_negAngle = @system(x' = Aext * x, x ∈ X)

    # deadzone
    A, b = get_dynamics(0., -α, u, n, θ)
    X = HPolyhedron([HalfSpace(sparsevec([1], [-1.], n+1), α),  # x1 >= -α
                     HalfSpace(sparsevec([1], [1.], n+1), α)])  # x1 <= α
    Aext = add_dimension(A, 1)
    Aext[1:n, n+1] .= b
    m_deadzone = @system(x' = Aext * x, x ∈ X)

    # posAngle
    A, b = get_dynamics(kₛ, α, u, n, θ)
    X = HalfSpace(sparsevec([1], [-1.], n+1), -α)  # x1 >= α
    Aext = add_dimension(A, 1)
    Aext[1:n, n+1] .= b
    m_posAngle = @system(x' = Aext * x, x ∈ X)

    # negAngleInit
    A, b = get_dynamics(kₛ, -α, -u, n, θ)
    X = HalfSpace(sparsevec([n], [1.], n+1), t_init)  # t <= t_init
    Aext = add_dimension(A, 1)
    Aext[1:n, n+1] .= b
    m_negAngleInit = @system(x' = Aext * x, x ∈ X)

    # modes
    modes = [m_negAngle, m_deadzone, m_posAngle, m_negAngleInit]

    # transition negAngleInit -> negAngle
    add_transition!(automaton, 4, 1, 1)
    guard = HalfSpace(sparsevec([n], [-1.], n+1), -t_init)  # t >= t_init
    r_41 = ConstrainedIdentityMap(n, guard)

    # transition negAngle -> deadzone
    add_transition!(automaton, 1, 2, 2)
    guard = HalfSpace(sparsevec([1], [-1.], n+1), α)  # x1 >= -α
    r_12 = ConstrainedIdentityMap(n, guard)

    # transition deadzone -> posAngle
    add_transition!(automaton, 2, 3, 3)
    guard = HalfSpace(sparsevec([1], [-1.], n+1), -α)  # x1 >= α
    r_23 = ConstrainedIdentityMap(n, guard)

    # transition annotations
    resetmaps = [r_41, r_12, r_23]

    # switching
    switchings = [HybridSystems.AutonomousSwitching()]

    H = HybridSystem(automaton, modes, resetmaps, switchings)

    # initial condition TODO check initial mode
    X0 = get_initial_condition(n, X0_scale)
    @assert isa(X0, Zonotope)
    X0 = Zonotope(vcat(X0.center, 1.0), vcat(X0.generators, 0.0))
    initial_condition = [(4, X0)]

    return InitialValueProblem(H, initial_condition)
end

const IA = ReachabilityAnalysis.IntervalArithmetic
using ReachabilityAnalysis: post

# requirement
function req(sol)
    idx = findall(x -> 2.0 ∈ tspan(x), sol.Fk.u)
    !isempty(idx) && all(x -> location(x) == 3, sol.Fk.u[idx])
end

function _solve_powertrain(prob)

    # solve mode 4
    p1 = IVP(mode(prob.s, 4), prob.x0[1][2])
    alg = GLGM06(δ=0.001, approx_model=Forward(setops=:zono))
    sol1F = post(alg, p1, IA.Interval(0, 2.0))
    sol1F.ext[:loc_id] = 4
    aux = sol1F(t_init) # clustered sets
    aux2 = overapproximate(ConvexHullArray([set(Xi) for Xi in aux]), Zonotope)

    # solve mode 1
    p2 = IVP(mode(prob.s, 1), aux2)
    Δt0 = tspan(aux)
    tprev = tstart(aux)
    sol2F = post(alg, p2, IA.Interval(0, 2.0 - tprev); Δt0=Δt0)
    sol2F.ext[:loc_id] = 1

    auxvec = zeros(dim(sol2F)) # = n + 1
    auxvec[1] = 1.0
    auxh = Hyperplane(auxvec, -α) # x₁ = -0.03
    jump_rset_idx = findall(R -> !isdisjoint(R, auxh), sol2F)
    sol2F_rset = view(sol2F, jump_rset_idx)
    aux3 = overapproximate(ConvexHullArray([set(Xi) for Xi in sol2F_rset]), Zonotope)

    # solve mode 2
    p3 = IVP(mode(prob.s, 2), aux3)
    alg2 = GLGM06(δ=0.001, approx_model=Forward(setops=:zono), disjointness_method=BoxEnclosure())
    Δt0 = tspan(sol2F_rset)
    tprev = tstart(sol2F_rset)
    sol3F = post(alg2, p3, IA.Interval(0, 2.0 - tprev); Δt0=Δt0)
    sol3F.ext[:loc_id] = 2
    auxh = Hyperplane(auxvec, α) # x₁ = 0.03
    jump_rset_idx = findall(R -> !isdisjoint(R, auxh), sol3F)
    sol3F_rset = view(sol3F, jump_rset_idx)
    aux4 = overapproximate(ConvexHullArray([set(Xi) for Xi in sol3F_rset]), Zonotope)

    # solve mode 3
    p4 = [IVP(mode(prob.s, 3), set(Xi)) for Xi in sol3F_rset]
    Δt0 = tspan(sol3F_rset)
    tprev = tstart(sol3F_rset)
    sol4F = [post(alg, p4i, IA.Interval(0, 2.0 - tprev); Δt0=Δt0) for p4i in p4]
    for F in sol4F
        F.ext[:loc_id] = 3
    end

    return HybridFlowpipe(reduce(vcat, (sol1F, sol2F, sol3F, sol4F)))
end
