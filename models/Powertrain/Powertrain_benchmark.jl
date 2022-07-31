using BenchmarkTools, Plots, Plots.PlotMeasures, LaTeXStrings
import LazySets
using BenchmarkTools: minimum, median

SUITE = BenchmarkGroup()
model = "POWERTRAIN"
cases = ["DTN01", "DTN02", "DTN03", "DTN04", "DTN05", "DTN06"]

SUITE[model] = BenchmarkGroup()

include("Powertrain.jl")
validation = []

LazySets.deactivate_assertions()

# ----------------------------------------
#  DTN01 (dense time)
# ----------------------------------------
prob_DTN01 = powertrain_homog(θ=1)
@time sol_DTN01 = _solve_powertrain(prob_DTN01) # warm-up run
property = req(sol_DTN01)
push!(validation, Int(property))
SUITE[model][cases[1]] = @benchmarkable _solve_powertrain($prob_DTN01)

# ----------------------------------------
#  DTN02 (dense time)
# ----------------------------------------
prob_DTN02 = powertrain_homog(θ=2)
@time sol_DTN02 = _solve_powertrain(prob_DTN02) # warm-up run
property = req(sol_DTN02)
push!(validation, Int(property))
SUITE[model][cases[2]] = @benchmarkable _solve_powertrain($prob_DTN02)

# ----------------------------------------
#  DTN03 (dense time)
# ----------------------------------------
prob_DTN03 = powertrain_homog(θ=3)
@time sol_DTN03 = _solve_powertrain(prob_DTN03) # warm-up run
property = req(sol_DTN03)
push!(validation, Int(property))
SUITE[model][cases[3]] = @benchmarkable _solve_powertrain($prob_DTN03)

# ----------------------------------------
#  DTN04 (dense time)
# ----------------------------------------
prob_DTN04 = powertrain_homog(θ=4)
@time sol_DTN04 = _solve_powertrain(prob_DTN04) # warm-up run
property = req(sol_DTN04)
push!(validation, Int(property))
SUITE[model][cases[4]] = @benchmarkable _solve_powertrain($prob_DTN04)

# ----------------------------------------
#  DTN05 (dense time)
# ----------------------------------------
prob_DTN05 = powertrain_homog(θ=5)
@time sol_DTN05 = _solve_powertrain(prob_DTN05) # warm-up run
property = req(sol_DTN05)
push!(validation, Int(property))
SUITE[model][cases[5]] = @benchmarkable _solve_powertrain($prob_DTN05)

# ----------------------------------------
#  DTN06 (dense time)
# ----------------------------------------
prob_DTN06 = powertrain_homog(θ=6)
@time sol_DTN06 = _solve_powertrain(prob_DTN06) # warm-up run
property = req(sol_DTN06)
push!(validation, Int(property))
SUITE[model][cases[6]] = @benchmarkable _solve_powertrain($prob_DTN06)

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
   print(io, "JuliaReach, $model, $c, $(validation[i]), $(runtimes[c])\n")
end

# ==============================================================================
# Plot
# ==============================================================================

fig = Plots.plot()
Plots.plot!(fig, sol_DTN03, vars=(1, 3),
           color=:blue, alpha=0.5, lw=1.0, linecolor=:black,
           tickfont=font(30, "Times"), guidefontsize=45,
           xlab=L"x_1",
           ylab=L"x_3",
           xtick=([-0.1, 0.0, 0.1, 0.2],
                  [L"-0.1", L"0.0", L"-0.1", L"0.2"]),
           ytick=([20, 40, 60, 80],
                  [L"20", L"40", L"60", L"80"]),
           xlims=(-0.1, 0.2), ylims=(-5.0, 85),
           bottom_margin=-1mm, left_margin=-2mm, right_margin=6mm, top_margin=0mm,
           size=(1000, 1000))
Plots.vline!(fig, [-α], lc=:red, ls=:dash, lw=2, lab="")
Plots.vline!(fig, [α], lc=:red, ls=:dash, lw=2, lab="")
savefig(fig, joinpath(@__DIR__, "ARCH-COMP22-JuliaReach-Powertrain_DTN03"))

fig = Plots.plot()
Plots.plot!(fig, sol_DTN05, vars=(1, 3),
           color=:blue, alpha=0.5, lw=1.0, linecolor=:black,
           tickfont=font(30, "Times"), guidefontsize=45,
           xlab=L"x_1",
           ylab=L"x_3",
           xtick=([-0.1, 0.0, 0.1, 0.2],
                  [L"-0.1", L"0.0", L"-0.1", L"0.2"]),
           ytick=([20, 40, 60, 80],
                  [L"20", L"40", L"60", L"80"]),
           xlims=(-0.1, 0.2), ylims=(-5.0, 85),
           bottom_margin=-1mm, left_margin=-2mm, right_margin=6mm, top_margin=0mm,
           size=(1000, 1000))
Plots.vline!(fig, [-α], lc=:red, ls=:dash, lw=2, lab="")
Plots.vline!(fig, [α], lc=:red, ls=:dash, lw=2, lab="")
savefig(fig, joinpath(@__DIR__, "ARCH-COMP22-JuliaReach-Powertrain_DTN05"))

sol_DTN01 = nothing
sol_DTN02 = nothing
sol_DTN03 = nothing
sol_DTN04 = nothing
sol_DTN05 = nothing
sol_DTN06 = nothing

GC.gc()
