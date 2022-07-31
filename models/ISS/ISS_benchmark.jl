using BenchmarkTools, Plots, Plots.PlotMeasures, LaTeXStrings
import LazySets
using BenchmarkTools: minimum, median

SUITE = BenchmarkGroup()
model = "ISS"
cases = ["ISSF01-ISS01", "ISSF01-ISU01",
         "ISSF01-ISS01-discrete", "ISSF01-ISU01-discrete",
         "ISSC01-ISS02", "ISSC01-ISU02",
         "ISSC01-ISS02-discrete", "ISSC01-ISU02-discrete"]

SUITE[model] = BenchmarkGroup()

include("ISS.jl")
validation = []

LazySets.deactivate_assertions()

# ----------------------------------------
#  ISSF01 (dense time)
# ----------------------------------------
prob_ISSF01 = ISSF01()

dirs = CustomDirections([C3, -C3])

# ISSF01-ISS01
sol_ISSF01 = solve(prob_ISSF01, T=20.0, alg=LGG09(δ=6e-4, template=dirs, sparse=true, cache=false))
property = (-ρ(-C3, sol_ISSF01) >= -0.0007) && (ρ(C3, sol_ISSF01) <= 0.0007)
push!(validation, Int(property))
SUITE[model][cases[1]] = @benchmarkable solve($prob_ISSF01, T=20.0, alg=LGG09(δ=6e-4, template=$dirs, sparse=true, cache=false))

# ISSF01-ISU01
sol_ISSF01 = solve(prob_ISSF01, T=20.0, alg=LGG09(δ=6e-4, template=dirs, sparse=true, cache=false))
property = (-ρ(-C3, sol_ISSF01) >= -0.0005) && (ρ(C3, sol_ISSF01) <= 0.0005)
push!(validation, Int(property))
SUITE[model][cases[2]] = @benchmarkable solve($prob_ISSF01, T=20.0, alg=LGG09(δ=6e-4, template=$dirs, sparse=true, cache=false))


# ----------------------------------------
#  ISSF01 (discrete time)
# ----------------------------------------

dirs = CustomDirections([C3, -C3])

# ISSF01-ISS01
sol_ISSF01 = solve(prob_ISSF01, T=20.0, alg=LGG09(δ=0.01, template=dirs, sparse=true, cache=false, approx_model=NoBloating()))
property = (-ρ(-C3, sol_ISSF01) >= -0.0007) && (ρ(C3, sol_ISSF01) <= 0.0007)
push!(validation, Int(property))
SUITE[model][cases[3]] = @benchmarkable solve($prob_ISSF01, T=20.0, alg=LGG09(δ=6e-4, template=$dirs, sparse=true, cache=false))

# ISSF01-ISU01
sol_ISSF01 = solve(prob_ISSF01, T=20.0, alg=LGG09(δ=0.01, template=dirs, sparse=true, cache=false, approx_model=NoBloating()))
property = (-ρ(-C3, sol_ISSF01) >= -0.0005) && (ρ(C3, sol_ISSF01) <= 0.0005)
push!(validation, Int(property))
SUITE[model][cases[4]] = @benchmarkable solve($prob_ISSF01, T=20.0, alg=LGG09(δ=6e-4, template=$dirs, sparse=true, cache=false))


# ----------------------------------------
#  ISSC01 (dense time)
# ----------------------------------------
prob_ISSC01 = ISSC01()

dirs = CustomDirections([C3_ext, -C3_ext])

# ISSC01-ISS02
sol_ISSC01 = solve(prob_ISSC01, T=20.0, alg=LGG09(δ=0.01, template=dirs, sparse=true, cache=true))
property = (-ρ(-C3_ext, sol_ISSC01) >= -0.0005) && (ρ(C3_ext, sol_ISSC01) <= 0.0005)
push!(validation, Int(property))
SUITE[model][cases[5]] = @benchmarkable solve($prob_ISSC01, T=20.0, alg=LGG09(δ=0.01, template=$dirs, sparse=true, cache=true))

# ISSC01-ISU02
sol_ISSC01 = solve(prob_ISSC01, T=20.0, alg=LGG09(δ=0.01, template=dirs, sparse=true, cache=true))
property = (-ρ(-C3_ext, sol_ISSC01) >= -0.00017) && (ρ(C3_ext, sol_ISSC01) <= 0.00017)
push!(validation, Int(property))
SUITE[model][cases[6]] = @benchmarkable solve($prob_ISSC01, T=20.0, alg=LGG09(δ=0.01, template=$dirs, sparse=true, cache=true))


# ----------------------------------------
#  ISSC01 (discrete time)
# ----------------------------------------

prob_ISSC01 = ISSC01()

dirs = CustomDirections([C3_ext, -C3_ext])

# ISSC01-ISS02
sol_ISSC01 = solve(prob_ISSC01, T=20.0, alg=LGG09(δ=0.01, template=dirs, sparse=true, cache=true, approx_model=NoBloating()))
property = (-ρ(-C3_ext, sol_ISSC01) >= -0.0005) && (ρ(C3_ext, sol_ISSC01) <= 0.0005)
push!(validation, Int(property))
SUITE[model][cases[7]] = @benchmarkable solve($prob_ISSC01, T=20.0, alg=LGG09(δ=0.01, template=$dirs, sparse=true, cache=true))

# ISSC01-ISU02
sol_ISSC01 = solve(prob_ISSC01, T=20.0, alg=LGG09(δ=0.01, template=dirs, sparse=true, cache=true, approx_model=NoBloating()))
property = (-ρ(-C3_ext, sol_ISSC01) >= -0.00017) && (ρ(C3_ext, sol_ISSC01) <= 0.00017)
push!(validation, Int(property))
SUITE[model][cases[8]] = @benchmarkable solve($prob_ISSC01, T=20.0, alg=LGG09(δ=0.01, template=$dirs, sparse=true, cache=true))

sol_ISSF01 = nothing
sol_ISSC01 = nothing
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
    t = median(results[model][c]).time * 1e-9
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

dirs = CustomDirections([C3, -C3])
prob_ISSF01 = ISSF01()
sol_ISSF01 = solve(prob_ISSF01, T=20.0, alg=LGG09(δ=6e-4, template=dirs, sparse=true, cache=false))
πsol_ISSF01 = flatten(sol_ISSF01)

fig = Plots.plot()
Plots.plot!(fig, πsol_ISSF01, vars=(0, 1),
            linecolor=:blue, color=:blue, alpha=0.8,
            tickfont=font(30, "Times"), guidefontsize=45,
            xlab=L"t",
            ylab=L"y_{3}",
            xtick=[0, 5, 10, 15, 20.], ytick=[-0.0005, 0, 0.0005],
            xlims=(0., 20.), ylims=(-0.00075, 0.00075),
            yformatter = y -> string(floor(Int, y/1e-4)), annotate=(0.5, 7.5e-4, text("⋅10⁻⁴", :left, 30)),
            bottom_margin=-6mm, left_margin=-1mm, right_margin=4mm, top_margin=3mm,
            size=(1000, 1000))
Plots.hline!(fig, [-0.0005], lc=:red, ls=:dash, lw=2, lab="")
Plots.hline!(fig, [0.0005], lc=:red, ls=:dash, lw=2, lab="")
Plots.hline!(fig, [-0.0007], lc=:red, ls=:dash, lw=2, lab="")
Plots.hline!(fig, [0.0007], lc=:red, ls=:dash, lw=2, lab="")

savefig(fig, joinpath(@__DIR__, "ARCH-COMP22-JuliaReach-ISS-ISSF01"))

# -------------

dirs = CustomDirections([C3_ext, -C3_ext])
prob_ISSC01 = ISSC01()
sol_ISSC01 = solve(prob_ISSC01, T=20.0, alg=LGG09(δ=0.01, template=dirs, sparse=true, cache=false))
πsol_ISSC01 = flatten(sol_ISSC01)

fig = Plots.plot()
Plots.plot!(fig, πsol_ISSC01, vars=(0, 1),
            linecolor=:blue, color=:blue, alpha=0.8,
            tickfont=font(30, "Times"), guidefontsize=45,
            xlab=L"t",
            ylab=L"y_{3}",
            xtick=[0, 5, 10, 15, 20.], ytick=[-0.0002, -0.0001, 0.0, 0.0001, 0.0002],
            xlims=(0., 20.), ylims=(-0.0002, 0.0002),
            yformatter = y -> string(floor(Int, y/1e-4)), annotate=(0.5, 2e-4, text("⋅10⁻⁴", :left, 30)),
            bottom_margin=-6mm, left_margin=-1mm, right_margin=4mm, top_margin=3mm,
            size=(1000, 1000))
Plots.hline!(fig, [-0.00017], lc=:red, ls=:dash, lw=2, lab="")
Plots.hline!(fig, [0.00017], lc=:red, ls=:dash, lw=2, lab="")

savefig(fig, joinpath(@__DIR__, "ARCH-COMP22-JuliaReach-ISS-ISSC01"))
