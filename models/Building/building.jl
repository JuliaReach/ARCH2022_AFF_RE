using ReachabilityAnalysis, SparseArrays, JLD2

LazySets.set_ztol(Float64, 1e-14)

const x25 = [zeros(24); 1.0; zeros(23)]
const x25e = vcat(x25, 0.0)

building_path = joinpath(@__DIR__, "building.jld2")

function building_BLDF01()
    @load building_path A B
    n = size(A, 1)
    U = Interval(0.8, 1.0)
    S = @system(x' = Ax + Bu, u ∈ U, x ∈ Universe(n))

    # initial states
    center_X0 = [fill(0.000225, 10); fill(0.0, 38)]
    radius_X0 = [fill(0.000025, 10); fill(0.0, 14); 0.0001; fill(0.0, 23)]
    X0 = Hyperrectangle(center_X0, radius_X0)

    prob_BLDF01 = InitialValueProblem(S, X0)
end

using ReachabilityAnalysis: add_dimension

function building_BLDC01()
    @load building_path A B
    n = size(A, 1)
    U = Interval(0.8, 1.0)

    # initial states
    center_X0 = [fill(0.000225, 10); fill(0.0, 38)]
    radius_X0 = [fill(0.000025, 10); fill(0.0, 14); 0.0001; fill(0.0, 23)]
    X0 = Hyperrectangle(center_X0, radius_X0)

    Ae = add_dimension(A)
    Ae[1:n, end] = B
    prob_BLDC01 = @ivp(x' = Ae * x, x(0) ∈ X0 × U)
end
