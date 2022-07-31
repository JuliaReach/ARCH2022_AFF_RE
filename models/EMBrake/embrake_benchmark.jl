using BenchmarkTools, Plots, Plots.PlotMeasures, LaTeXStrings
import LazySets
using BenchmarkTools: minimum, median

SUITE = BenchmarkGroup()
model = "EMBRAKE"
cases = ["BRKDC01", "BRKDC01-D", "BRKNC01", "BRKNC01-D", "BRKNP01", "BRKNP01-D"]
SUITE[model] = BenchmarkGroup()

include("embrake.jl")
validation = []
times = []

LazySets.deactivate_assertions()

# ----------------------------------------
#  BRKDC01
# ----------------------------------------

prob_no_pv_no_jit = embrake_no_pv(ζ=0., Tsample=1e-4)
alg = GLGM06(δ=1e-7, max_order=1, static=true, dim=4, ngens=4, approx_model=Forward())
sol_no_pv_no_jit = solve(prob_no_pv_no_jit, max_jumps=1001, alg=alg)

# verify that specification holds
property = ρ(eₓ, sol_no_pv_no_jit) < x0
push!(validation, Int(property))

t = max_time(sol_no_pv_no_jit)
push!(times, trunc(t, digits=4))
println("maximum time that x < x0 , case $(cases[1]) : $t")

# benchmark
SUITE[model][cases[1]] = @benchmarkable solve($prob_no_pv_no_jit, max_jumps=1001, alg=$alg)

# ----------------------------------------
#  BRKDC01 (discrete-time)
# ----------------------------------------

alg = GLGM06(δ=1e-7, max_order=1, static=true, dim=4, ngens=4, approx_model=NoBloating())
sol_no_pv_no_jit_discrete = solve(prob_no_pv_no_jit, max_jumps=1001, alg=alg)

# verify that specification holds
property = ρ(eₓ, sol_no_pv_no_jit_discrete) < x0
push!(validation, Int(property))

t = max_time(sol_no_pv_no_jit_discrete)
push!(times, trunc(t, digits=4))
println("maximum time that x < x0 , case $(cases[2]) : $t")

# benchmark
SUITE[model][cases[2]] = @benchmarkable solve($prob_no_pv_no_jit, max_jumps=1001, alg=$alg)

sol_no_pv_no_jit_discrete = nothing
GC.gc()

# ----------------------------------------
#  BRKNC01
# ----------------------------------------

prob_no_pv_with_jit = embrake_no_pv(ζ=[-1e-8, 1e-7], Tsample=1e-4)
alg = GLGM06(δ=1e-8, max_order=1, static=true, dim=4, ngens=4, approx_model=Forward())
sol_no_pv_with_jit = solve(prob_no_pv_with_jit, max_jumps=1001, alg=alg)

# verify that specification holds
property = ρ(eₓ, sol_no_pv_with_jit) < x0
push!(validation, Int(property))

GC.gc()
t = max_time(sol_no_pv_with_jit)
push!(times, trunc(t, digits=4))
println("maximum time that x < x0 , case $(cases[3]) : $t")

# benchmark
SUITE[model][cases[3]] = @benchmarkable solve($prob_no_pv_with_jit, max_jumps=1001, alg=$alg)

sol_no_pv_with_jit = nothing
GC.gc()

# ----------------------------------------
#  BRKNC01 (discrete-time)
# ----------------------------------------

alg = GLGM06(δ=1e-8, max_order=1, static=true, dim=4, ngens=4, approx_model=NoBloating())
sol_no_pv_no_jit_discrete = solve(prob_no_pv_no_jit, max_jumps=1001, alg=alg)

# verify that specification holds
property = ρ(eₓ, sol_no_pv_no_jit_discrete) < x0
push!(validation, Int(property))

GC.gc()
t = max_time(sol_no_pv_no_jit_discrete)
push!(times, trunc(t, digits=4))
println("maximum time that x < x0 , case $(cases[4]) : $t")

# benchmark
SUITE[model][cases[4]] = @benchmarkable solve($prob_no_pv_no_jit, max_jumps=1001, alg=$alg)

sol_no_pv_no_jit_discrete = nothing
GC.gc()

# ----------------------------------------
#  BRKNP01
# ----------------------------------------

prob_pv_1_with_jit = embrake_pv_1(ζ=[-1e-8, 1e-7], Tsample=1e-4)
alg = ASB07(δ=1e-8, max_order=1, static=true, dim=4, ngens=4)
sol_pv_1_with_jit = solve(prob_pv_1_with_jit, max_jumps=1001, alg=alg)

# verify that specification holds
property = ρ(eₓ, sol_pv_1_with_jit) < x0
push!(validation, Int(property))

GC.gc()
t = max_time(sol_pv_1_with_jit)
push!(times, trunc(t, digits=4))
println("maximum time that x < x0 , case $(cases[5]) : $t")

# benchmark
SUITE[model][cases[5]] = @benchmarkable solve($prob_pv_1_with_jit, max_jumps=1001, alg=$alg)

sol_pv_1_with_jit = nothing
GC.gc()

# ----------------------------------------
#  BRKNP01 (discrete-time)
# ----------------------------------------

alg = ASB07(δ=1e-8, max_order=1, static=true, dim=4, ngens=4, approx_model=NoBloating(exp=:interval))
sol_pv_1_with_jit_discrete = solve(prob_pv_1_with_jit, max_jumps=1001, alg=alg)

# verify that specification holds
property = ρ(eₓ, sol_pv_1_with_jit_discrete) < x0
push!(validation, Int(property))

GC.gc()
t = max_time(sol_pv_1_with_jit_discrete)
push!(times, trunc(t, digits=4))
println("maximum time that x < x0 , case $(cases[6]) : $t")

# benchmark
SUITE[model][cases[6]] = @benchmarkable solve($prob_pv_1_with_jit, max_jumps=1001, alg=$alg)

sol_pv_1_with_jit_discrete = nothing
GC.gc()

# ==============================================================================
# Execute benchmarks and save benchmark results
# ==============================================================================

# tune parameters
tune!(SUITE)

# run the benchmarks
results = run(SUITE, verbose=true)

# return the sample with the smallest time value in each test
println("minimum time for each benchmark:\n", minimum(results))

# return the median for each test
println("median time for each benchmark:\n", median(results))

# export runtimes
runtimes = Dict()
for (i, c) in enumerate(cases)
    local t = median(results[model][c]).time * 1e-9
    runtimes[c] = round(t, digits=4)
end

if !@isdefined io
    io = stdout
end

for (i, c) in enumerate(cases)
    print(io, "JuliaReach, $model, $c, $(validation[i]), $(runtimes[c]), $(times[i])\n")
end

# ==============================================================================
# Plot
# ==============================================================================

polys = Vector{VPolygon{Float64, Vector{Float64}}}()
for fp in sol_no_pv_no_jit
    for (j, X) in enumerate(fp)
        sfpos, sfneg = ρ(eₓ, X), ρ(-eₓ, X)
        dt = tspan(fp, j)
        ti, tf = inf(dt), sup(dt)
        p = VPolygon([[ti, sfpos], [ti, -sfneg], [tf, sfpos], [tf, -sfneg]])
        push!(polys, p)
    end
end

fig = Plots.plot()

plot!(fig, [X for X in polys[1:500:end]], color=:black, lw=1.0, linecolor=:black,
        tickfont=font(30, "Times"), guidefontsize=45,
        xlab=L"t",
        ylab=L"x",
        xtick=[0.025, 0.05, 0.075, 0.1], ytick=[0.0, 0.01, 0.02, 0.03, 0.04, 0.05],
        xlims=(0.0, 0.1), ylims=(0.0, 0.05),
        bottom_margin=-6mm, left_margin=-3mm, right_margin=12mm, top_margin=3mm,
        size=(1000, 1000))
Plots.hline!(fig, [x0], lc=:red, ls=:dash, lw=2, lab="")

savefig(fig, joinpath(@__DIR__, "ARCH-COMP22-JuliaReach-EMBrake.png"))
# savefig(fig, joinpath(@__DIR__, "ARCH-COMP22-JuliaReach-EMBrake.pdf"))

sol_no_pv_no_jit = nothing
GC.gc()
