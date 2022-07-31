#=
This model is taken from [1]. The theoretical derivations of the analytic solution can be found in [2].

[1] Malakiyeh, Mohammad Mahdi, Saeed Shojaee, and Klaus-Jürgen Bathe. "The Bathe time integration method revisited for prescribing desired numerical dissipation." Computers & Structures 212 (2019): 289-298.

[2] Mechanical Vibrations, Gerardin et al, page 250-251.
=#

using ReachabilityAnalysis, LinearAlgebra, LazySets
using ReachabilityAnalysis: solve, discretize
using LazySets.Arrays
using SparseArrays
using SpaceExParser, MathematicalSystems

const CLCCS = ConstrainedLinearControlContinuousSystem

include("clamped_parse.jl")

LazySets.set_ztol(Float64, 1e-15)
LazySets.set_atol(Float64, 1e-15)
LazySets.set_rtol(Float64, 1e-15)

function clamped(model::String; mode=nothing, N=nothing, Z=nothing)
    # Z is the inputs assumption, either "C" or "F"
    if isnothing(mode) || isnothing(N) || isnothing(Z)
        mode, N, Z = parse_clamped_beam(model)
    end
    A = sparse(mode.A)
    ΔF0 = Interval(0.99, 1.01)

    if Z == "C"
        S = LinearContinuousSystem(A)
        X0 = Singleton(zeros(2N)) × ΔF0
    elseif Z == "F"
        S = CLCCS(A, mode.B, Universe(2N), ΔF0)
        X0 = Singleton(zeros(2N))
    end
    prob = InitialValueProblem(S, X0)
end
