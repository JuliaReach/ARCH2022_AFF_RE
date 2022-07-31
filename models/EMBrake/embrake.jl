using ReachabilityAnalysis, SparseArrays

LazySets.set_rtol(Float64, 1e-13)
LazySets.set_ztol(Float64, 1e-13)

const eₓ = sparsevec([2], [1.], 4)
const x0 = 0.05

# model without parameter variation
function embrake_no_pv(; Tsample=1.E-4, ζ=1e-6, x0=0.05)
    # model's constants
    L = 1.e-3
    KP = 10000.
    KI = 1000.
    R = 0.5
    K = 0.02
    drot = 0.1
    i = 113.1167

    # state variables: [I, x, xe, xc]
    A = Matrix([-(R+K^2/drot)/L 0      KP/L  KI/L;
                          K/i/drot      0      0      0;
                          0             0      0      0;
                          0             0      0      0])

    EMbrake = @system(x' = Ax)

    # initial conditions
    I₀  = Singleton([0.0])
    x₀  = Singleton([0.0])
    xe₀ = Singleton([0.0])
    xc₀ = Singleton([0.0])
    X₀ = I₀ × x₀ × xe₀ × xc₀

    # reset map
    Ar = sparse([1, 2, 3, 4, 4], [1, 2, 2, 2, 4], [1., 1., -1., -Tsample, 1.], 4, 4)
    br = sparsevec([3, 4], [x0, Tsample*x0], 4)
    reset_map(X) = Ar * X + br

    # hybrid system with clocked linear dynamics
    ha = HACLD1(EMbrake, reset_map, Tsample, ζ)
    return IVP(ha, X₀)
end

# model with parameter variation changing only 1 coefficient, corresponds to the Flow* settings in [SO15]
function embrake_pv_1(; Tsample=1.E-4, ζ=1e-6, Δ=3.0, x0=0.05)
    # model's constants
    L = 1.e-3
    KP = 10000.
    KI = 1000.
    R = 0.5
    K = 0.02
    drot = 0.1
    i = 113.1167
    p = 504. + (-Δ .. Δ)

    # state variables: [I, x, xe, xc]
    A = IntervalMatrix([-p            0      KP/L   KI/L;
                        K/i/drot      0      0      0;
                        0             0      0      0;
                        0             0      0      0])

    EMbrake = @system(x' = Ax)

    # initial conditions
    I₀  = Singleton([0.0])
    x₀  = Singleton([0.0])
    xe₀ = Singleton([0.0])
    xc₀ = Singleton([0.0])
    X₀ = I₀ × x₀ × xe₀ × xc₀

    # reset map
    Ar = sparse([1, 2, 3, 4, 4], [1, 2, 2, 2, 4], [1., 1., -1., -Tsample, 1.], 4, 4)
    br = sparsevec([3, 4], [x0, Tsample*x0], 4)
    reset_map(X) = Ar * X + br

    # hybrid system with clocked linear dynamics
    ha = HACLD1(EMbrake, reset_map, Tsample, ζ)
    return IVP(ha, X₀)
end

# maximum time by which the solution is below the x0 threshold
function max_time(sol)
    sfpos = []
    times = []

    for fp in sol
        for (j, X) in enumerate(fp)
            push!(sfpos, ρ(eₓ, X))
            push!(times, tspan(fp, j))
        end
    end
    idx = findfirst(x -> x > x0, sfpos)

    if isnothing(idx)
        return tend(sol)
    else
        return sup(times[idx-1])
    end
end
