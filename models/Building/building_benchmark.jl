using BenchmarkTools, Plots, Plots.PlotMeasures, LaTeXStrings
import LazySets
using BenchmarkTools: minimum, median

SUITE = BenchmarkGroup()
model = "BUILDING"
cases = ["BLDF01-BDS01", "BLDF01-BDU01", "BLDF01-BDU02",
         "BLDF01-BDS01-discrete", "BLDF01-BDU01-discrete", "BLDF01-BDU02-discrete",
         "BLDC01-BDS01", "BLDC01-BDU01", "BLDC01-BDU02",
         "BLDC01-BDS01-discrete", "BLDC01-BDU01-discrete", "BLDC01-BDU02-discrete"]

SUITE[model] = BenchmarkGroup()

include("building.jl")
validation = []

LazySets.deactivate_assertions()

# ----------------------------------------
#  BLDF01 (dense time)
# ----------------------------------------
prob_BLDF01 = building_BLDF01()

# BLDF01 - BDS01
sol_BLDF01 = solve(prob_BLDF01, T=20.0, alg=LGG09(δ=0.004, template=x25))
property = ρ(x25, sol_BLDF01) <= 0.0051
push!(validation, Int(property))
SUITE[model][cases[1]] = @benchmarkable solve($prob_BLDF01, T=20.0, alg=LGG09(δ=0.004, template=$x25))

# BLDF01 - BDU01
sol_BLDF01 = solve(prob_BLDF01, T=20.0, alg=LGG09(δ=0.004, template=x25))
property = ρ(x25, sol_BLDF01) <= 4e-3
push!(validation, Int(property))
SUITE[model][cases[2]] = @benchmarkable solve($prob_BLDF01, T=20.0, alg=LGG09(δ=0.004, template=$x25))

# BLDF01 - BDU02
sol_BLDF01 = solve(prob_BLDF01, T=20.0, alg=LGG09(δ=0.004, template=x25))
property = ρ(x25, sol_BLDF01(20.0)) <= -0.78e-3
push!(validation, Int(property))
SUITE[model][cases[3]] = @benchmarkable solve($prob_BLDF01, T=20.0, alg=LGG09(δ=0.004, template=$x25))

# ----------------------------------------
#  BLDF01 (discrete time)
# ----------------------------------------

# BLDF01 - BDS01
sol_BLDF01 = solve(prob_BLDF01, T=20.0, alg=LGG09(δ=0.01, template=x25, approx_model=NoBloating()))
property = ρ(x25, sol_BLDF01) <= 0.0051
push!(validation, Int(property))
SUITE[model][cases[4]] = @benchmarkable solve($prob_BLDF01, T=20.0, alg=LGG09(δ=0.01, template=$x25, approx_model=NoBloating()))

# BLDF01 - BDU01
sol_BLDF01 = solve(prob_BLDF01, T=20.0, alg=LGG09(δ=0.01, template=x25, approx_model=NoBloating()))
property = ρ(x25, sol_BLDF01) <= 4e-3
push!(validation, Int(property))
SUITE[model][cases[5]] = @benchmarkable solve($prob_BLDF01, T=20.0, alg=LGG09(δ=0.01, template=$x25, approx_model=NoBloating()))

# BLDF01 - BDU02
sol_BLDF01 = solve(prob_BLDF01, T=20.0, alg=LGG09(δ=0.01, template=x25, approx_model=NoBloating()))
property = ρ(x25, sol_BLDF01(20.0)) <= -0.78e-3
push!(validation, Int(property))
SUITE[model][cases[6]] = @benchmarkable solve($prob_BLDF01, T=20.0, alg=LGG09(δ=0.01, template=$x25, approx_model=NoBloating()))

# ----------------------------------------
#  BLDC01 (dense time)
# ----------------------------------------
prob_BLDC01 = building_BLDC01()

# BLDC01 - BDS01
sol_BLDC01 = solve(prob_BLDC01, T=20.0, alg=LGG09(δ=0.005, template=x25e))
property = ρ(x25e, sol_BLDC01) <= 0.0051
push!(validation, Int(property))
SUITE[model][cases[7]] = @benchmarkable solve($prob_BLDC01, T=20.0, alg=LGG09(δ=0.005, template=$x25e))

# BLDC01 - BDU01
sol_BLDC01 = solve(prob_BLDC01, T=20.0, alg=LGG09(δ=0.005, template=x25e))
property = ρ(x25e, sol_BLDC01) <= 4e-3
push!(validation, Int(property))
SUITE[model][cases[8]] = @benchmarkable solve($prob_BLDC01, T=20.0, alg=LGG09(δ=0.005, template=$x25e))

# BLDC01 - BDU02
sol_BLDC01 = solve(prob_BLDC01, T=20.0, alg=LGG09(δ=0.005, template=x25e))
property = ρ(x25e, sol_BLDC01(20.0)) <= -0.78e-3
push!(validation, Int(property))
SUITE[model][cases[9]] = @benchmarkable solve($prob_BLDC01, T=20.0, alg=LGG09(δ=0.005, template=$x25e))

# ----------------------------------------
#  BLDC01 (discrete time)
# ----------------------------------------

# BLDC01 - BDS01
sol_BLDC01 = solve(prob_BLDC01, T=20.0, alg=LGG09(δ=0.01, template=x25e, approx_model=NoBloating()))
property = ρ(x25e, sol_BLDC01) <= 0.0051
push!(validation, Int(property))
SUITE[model][cases[10]] = @benchmarkable solve($prob_BLDC01, T=20.0, alg=LGG09(δ=0.01, template=$x25e, approx_model=NoBloating()))

# BLDC01 - BDU01
sol_BLDC01 = solve(prob_BLDC01, T=20.0, alg=LGG09(δ=0.01, template=x25e, approx_model=NoBloating()))
property = ρ(x25e, sol_BLDC01) <= 4e-3
push!(validation, Int(property))
SUITE[model][cases[11]] = @benchmarkable solve($prob_BLDC01, T=20.0, alg=LGG09(δ=0.01, template=$x25e, approx_model=NoBloating()))

# BLDC01 - BDU02
sol_BLDC01 = solve(prob_BLDC01, T=20.0, alg=LGG09(δ=0.01, template=x25e, approx_model=NoBloating()))
property = ρ(x25e, sol_BLDC01(20.0)) <= -0.78e-3
push!(validation, Int(property))
SUITE[model][cases[12]] = @benchmarkable solve($prob_BLDC01, T=20.0, alg=LGG09(δ=0.01, template=$x25e, approx_model=NoBloating()))


sol_BLDF01 = nothing
sol_BLDC01 = nothing
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
    print(io, "JuliaReach, $model, $c, $(validation[i]), $(runtimes[c])\n")
end

# ==============================================================================
# Plot
# ==============================================================================
#
dirs = CustomDirections([x25, -x25])

sol_BLDF01 = solve(prob_BLDF01, T=20.0, alg=LGG09(δ=0.004, template=dirs))
fig = Plots.plot()
Plots.plot!(fig, sol_BLDF01, vars=(0, 25),
            color=:blue, alpha=1.0, lw=1.0, linecolor=:blue,
            tickfont=font(30, "Times"), guidefontsize=45,
            xlab=L"t",
            ylab=L"x",
            xtick=[0, 5, 10, 15, 20], ytick=[-0.006, -0.003, 0.0, 0.003, 0.006],
            xlims=(0.0, 20.0), ylims=(-0.007, 0.006),
            yformatter = y -> string(floor(Int, y/1e-3)), annotate=(0.5, 6e-3, text("⋅10⁻³", :left, 30)),
            bottom_margin=0mm, left_margin=0mm, right_margin=4mm, top_margin=3mm,
            size=(1000, 1000))
Plots.hline!(fig, [0.0051], lc=:red, ls=:dash, lw=2, lab="")
Plots.hline!(fig, [4e-3], lc=:red, ls=:dash, lw=2, lab="")
Plots.hline!(fig, [-0.78e-3], lc=:red, ls=:dash, lw=2, lab="")
savefig(fig, joinpath(@__DIR__, "ARCH-COMP22-JuliaReach-Building_BLDF01_t20"))

sol_BLDF01 = solve(prob_BLDF01, T=1.0, alg=LGG09(δ=0.001, template=dirs))
fig = Plots.plot()
Plots.plot!(fig, sol_BLDF01, vars=(0, 25),
            color=:blue, alpha=1.0, lw=1.0, linecolor=:blue,
            tickfont=font(30, "Times"), guidefontsize=45,
            xlab=L"t",
            ylab=L"x",
            xtick=[0, 0.25, 0.50, 0.75, 1.0], ytick=[-0.006, -0.003, 0.0, 0.003, 0.006],
            xlims=(0.0, 1.0), ylims=(-0.007, 0.006),
            yformatter = y -> string(floor(Int, y/1e-3)), annotate=(0.025, 6e-3, text("⋅10⁻³", :left, 30)),
            bottom_margin=0mm, left_margin=0mm, right_margin=8mm, top_margin=3mm,
            size=(1000, 1000))
Plots.hline!(fig, [0.0051], lc=:red, ls=:dash, lw=2, lab="")
Plots.hline!(fig, [4e-3], lc=:red, ls=:dash, lw=2, lab="")
Plots.hline!(fig, [-0.78e-3], lc=:red, ls=:dash, lw=2, lab="")
savefig(fig, joinpath(@__DIR__, "ARCH-COMP22-JuliaReach-Building_BLDF01_t1"))

dirs = CustomDirections([x25e, -x25e])

sol_BLDC01 = solve(prob_BLDC01, T=20.0, alg=LGG09(δ=0.004, template=dirs))
fig = Plots.plot()
Plots.plot!(fig, sol_BLDC01, vars=(0, 25),
            color=:blue, alpha=1.0, lw=1.0, linecolor=:blue,
            tickfont=font(30, "Times"), guidefontsize=45,
            xlab=L"t",
            ylab=L"x",
            xtick=[0, 5, 10, 15, 20], ytick=[-0.006, -0.003, 0.0, 0.003, 0.006],
            xlims=(0.0, 20.0), ylims=(-0.007, 0.006),
            yformatter = y -> string(floor(Int, y/1e-3)), annotate=(0.5, 6e-3, text("⋅10⁻³", :left, 30)),
            bottom_margin=0mm, left_margin=0mm, right_margin=4mm, top_margin=3mm,
            size=(1000, 1000))
Plots.hline!(fig, [0.0051], lc=:red, ls=:dash, lw=2, lab="")
Plots.hline!(fig, [4e-3], lc=:red, ls=:dash, lw=2, lab="")
Plots.hline!(fig, [-0.78e-3], lc=:red, ls=:dash, lw=2, lab="")
savefig(fig, joinpath(@__DIR__, "ARCH-COMP22-JuliaReach-Building_BLDC01_t20"))

sol_BLDC01 = solve(prob_BLDC01, T=1.0, alg=LGG09(δ=0.001, template=dirs))
fig = Plots.plot()
Plots.plot!(fig, sol_BLDF01, vars=(0, 25),
            color=:blue, alpha=1.0, lw=1.0, linecolor=:blue,
            tickfont=font(30, "Times"), guidefontsize=45,
            xlab=L"t",
            ylab=L"x",
            xtick=[0, 0.25, 0.50, 0.75, 1.0], ytick=[-0.006, -0.003, 0.0, 0.003, 0.006],
            xlims=(0.0, 1.0), ylims=(-0.007, 0.006),
            yformatter = y -> string(floor(Int, y/1e-3)), annotate=(0.025, 6e-3, text("⋅10⁻³", :left, 30)),
            bottom_margin=0mm, left_margin=0mm, right_margin=8mm, top_margin=3mm,
            size=(1000, 1000))
Plots.hline!(fig, [0.0051], lc=:red, ls=:dash, lw=2, lab="")
Plots.hline!(fig, [4e-3], lc=:red, ls=:dash, lw=2, lab="")
Plots.hline!(fig, [-0.78e-3], lc=:red, ls=:dash, lw=2, lab="")
savefig(fig, joinpath(@__DIR__, "ARCH-COMP22-JuliaReach-Building_BLDC01_t1"));
