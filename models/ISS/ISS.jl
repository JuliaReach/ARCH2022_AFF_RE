using ReachabilityAnalysis, JLD2
using ReachabilityAnalysis: add_dimension

LazySets.set_ztol(Float64, 1e-15)
ISS_path = joinpath(@__DIR__, "ISS.jld2")

@load ISS_path C
const C3 = C[3, :] # variable y₃
const C3_ext = vcat(C3, fill(0.0, 3))

function ISSF01()
    @load ISS_path A B

    U = Hyperrectangle(low=[0.0, 0.8, 0.9], high=[0.1, 1., 1.]);  # input set
    X0 = BallInf(zeros(size(A, 1)), 0.0001)  # -0.0001 <= xi <= 0.0001 for all i
    prob_ISSF01 = @ivp(x' = A*x + B*u, x(0) ∈ X0, u ∈ U, x ∈ Universe(270))
end

function ISSC01()
    @load ISS_path A B

    Aext = add_dimension(A, 3)
    Aext[1:270, 271:273] = B

    U = Hyperrectangle(low=[0.0, 0.8, 0.9], high=[0.1, 1., 1.])  # input set
    X0 = BallInf(zeros(size(A, 1)), 0.0001)  # -0.0001 <= xi <= 0.0001 for all i
    X0 = X0 * U
    prob_ISSC01 = @ivp(x' = Aext*x, x(0) ∈ X0)
end
