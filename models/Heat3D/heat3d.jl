using ReachabilityAnalysis, SparseArrays, ExponentialUtilities, JLD2
using LazySets.Arrays
using ReachabilityAnalysis: reach_homog_krylov_LGG09!

# For xc we add +1 because Julia has 1based indexing in contrast to Python which is 0 based

# HEAT01 5x5x5 grid
function heat01(; δ=0.02) # step size
    @load joinpath(@__DIR__, "HEAT01.jld2") A xind

    n = size(A, 1)
    x0ind = xind.nzind
    xc = 62 + 1

    c = sparsevec(x0ind, fill(1.0, length(x0ind)), n)
    r = sparsevec(x0ind, fill(0.1, length(x0ind)), n)
    Ω₀ = Hyperrectangle(c, r)

    Aᵀδ = transpose(A .* δ)
    ℓ = SingleEntryVector(xc, n, 1.0)

    return A, Aᵀδ, Ω₀, ℓ
end

# HEAT02 10x10x10 grid
function heat02(; δ=0.02) # step size
    @load joinpath(@__DIR__, "HEAT02.jld2") A xind

    n = size(A, 1)
    x0ind = xind.nzind
    xc = 555 + 1

    c = sparsevec(x0ind, fill(1.0, length(x0ind)), n)
    r = sparsevec(x0ind, fill(0.1, length(x0ind)), n)
    Ω₀ = Hyperrectangle(c, r)

    Aᵀδ = transpose(A .* δ)
    ℓ = SingleEntryVector(xc, n, 1.0)

    return A, Aᵀδ, Ω₀, ℓ
end

# HEAT03 20x20x20 grid
function heat03(; δ=0.02) # step size
    @load joinpath(@__DIR__, "HEAT03.jld2") A xind

    n = size(A, 1)
    x0ind = xind.nzind
    xc = 4210 + 1

    c = sparsevec(x0ind, fill(1.0, length(x0ind)), n)
    r = sparsevec(x0ind, fill(0.1, length(x0ind)), n)
    Ω₀ = Hyperrectangle(c, r)

    Aᵀδ = transpose(A .* δ)
    ℓ = SingleEntryVector(xc, n, 1.0)

    return A, Aᵀδ, Ω₀, ℓ
end

# HEAT04 50x50x50 grid
function heat04(; δ=0.02) # step size
    @load joinpath(@__DIR__, "HEAT04.jld2") A xind

    n = size(A, 1)
    x0ind = xind.nzind
    xc = 63775 + 1

    c = sparsevec(x0ind, fill(1.0, length(x0ind)), n)
    r = sparsevec(x0ind, fill(0.1, length(x0ind)), n)
    Ω₀ = Hyperrectangle(c, r)

    Aᵀδ = transpose(A .* δ)
    ℓ = SingleEntryVector(xc, n, 1.0)

    return A, Aᵀδ, Ω₀, ℓ
end
