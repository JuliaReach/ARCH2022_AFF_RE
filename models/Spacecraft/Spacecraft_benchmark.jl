using BenchmarkTools, Plots, Plots.PlotMeasures, LaTeXStrings
import LazySets
using BenchmarkTools: minimum, median

SUITE = BenchmarkGroup()
model = "SPACECRAFT"
cases = ["NA01", "NA01-discrete",
         "A01", "A01-discrete",
         "A02", "A02-discrete",
         "A03", "A03-discrete",
         "A04", # "A04-discrete",
         #"A05", "A05-discrete",
         #"A06", "A06-discrete",
         #"A07", "A07-discrete",
         #"A08", "A08-discrete",
         "U01", # "U01-discrete",
         "U02", "U02-discrete"]

SUITE[model] = BenchmarkGroup()

include("Spacecraft.jl")
validation = []

LazySets.deactivate_assertions()

boxdirs = BoxDirections{Float64, Vector{Float64}}(5)

# ----------------------------------------
#  NA01 (dense time)
# ----------------------------------------

prob_NA01 = spacecraft(abort_time=-1.)

sol_NA01 = solve(prob_NA01, alg=BOX(δ=0.04),
                 clustering_method=LazyClustering(),
                 intersection_method=TemplateHullIntersection(boxdirs),
                 intersect_source_invariant=false,
                 tspan = (0.0 .. 300.0))
property = SR02_specification(sol_NA01)
push!(validation, Int(property))
SUITE[model][cases[1]] = @benchmarkable solve($prob_NA01, alg=BOX(δ=0.04),
                 clustering_method=LazyClustering(),
                 intersection_method=TemplateHullIntersection($boxdirs),
                 intersect_source_invariant=false,
                 tspan = (0.0 .. 300.0))

# ----------------------------------------
#  NA01 (discrete time)
# ----------------------------------------

sol_NA01 = solve(prob_NA01, alg=BOX(δ=0.1, approx_model=NoBloating()),
                 clustering_method=LazyClustering(),
                 intersection_method=TemplateHullIntersection(boxdirs),
                 intersect_source_invariant=false,
                 tspan = (0.0 .. 300.0))
property = SR02_specification(sol_NA01)
push!(validation, Int(property))
SUITE[model][cases[2]] = @benchmarkable solve($prob_NA01, alg=BOX(δ=0.1, approx_model=NoBloating()),
                 clustering_method=LazyClustering(),
                 intersection_method=TemplateHullIntersection($boxdirs),
                 intersect_source_invariant=false,
                 tspan = (0.0 .. 300.0))

 # ----------------------------------------
 #  A01 (dense time)
 # ----------------------------------------

 prob_A01 = spacecraft(abort_time=120.)

 sol_A01 = solve(prob_A01, alg=BOX(δ=0.04),
                 clustering_method=LazyClustering(3),
                 intersection_method=TemplateHullIntersection(boxdirs),
                 intersect_source_invariant=false,
                 tspan = (0.0 .. 300.0))
 property = SR02_specification(sol_A01)
 push!(validation, Int(property))
 SUITE[model][cases[3]] = @benchmarkable solve($prob_A01, alg=BOX(δ=0.04),
                  clustering_method=LazyClustering(3),
                  intersection_method=TemplateHullIntersection($boxdirs),
                  intersect_source_invariant=false,
                  tspan = (0.0 .. 300.0))

 # ----------------------------------------
 #  A01 (discrete time)
 # ----------------------------------------

 sol_A01 = solve(prob_A01, alg=BOX(δ=0.1, approx_model=NoBloating()),
                 clustering_method=LazyClustering(3),
                 intersection_method=TemplateHullIntersection(boxdirs),
                 intersect_source_invariant=false,
                 tspan = (0.0 .. 300.0))
 property = SR02_specification(sol_A01)
 push!(validation, Int(property))
 SUITE[model][cases[4]] = @benchmarkable solve($prob_A01, alg=BOX(δ=0.1, approx_model=NoBloating()),
                  clustering_method=LazyClustering(3),
                  intersection_method=TemplateHullIntersection($boxdirs),
                  intersect_source_invariant=false,
                  tspan = (0.0 .. 300.0))


# ----------------------------------------
#  A02 (dense time)
# ----------------------------------------

prob_A02 = spacecraft(abort_time=[120., 125.])

sol_A02 = solve(prob_A02, alg=BOX(δ=0.04),
              clustering_method=LazyClustering(3),
              intersection_method=TemplateHullIntersection(boxdirs),
              intersect_source_invariant=false,
              tspan = (0.0 .. 300.0))
property = SR02_specification(sol_A02)
push!(validation, Int(property))
SUITE[model][cases[5]] = @benchmarkable solve($prob_A02, alg=BOX(δ=0.04),
               clustering_method=LazyClustering(3),
               intersection_method=TemplateHullIntersection($boxdirs),
               intersect_source_invariant=false,
               tspan = (0.0 .. 300.0))

# ----------------------------------------
#  A02 (discrete time)
# ----------------------------------------

sol_A02 = solve(prob_A02, alg=BOX(δ=0.1, approx_model=NoBloating()),
              clustering_method=LazyClustering(3),
              intersection_method=TemplateHullIntersection(boxdirs),
              intersect_source_invariant=false,
              tspan = (0.0 .. 300.0))
property = SR02_specification(sol_A02)
push!(validation, Int(property))
SUITE[model][cases[6]] = @benchmarkable solve($prob_A02, alg=BOX(δ=0.1, approx_model=NoBloating()),
               clustering_method=LazyClustering(3),
               intersection_method=TemplateHullIntersection($boxdirs),
               intersect_source_invariant=false,
               tspan = (0.0 .. 300.0))


# ----------------------------------------
#  A03 (dense time)
# ----------------------------------------

prob_A03 = spacecraft(abort_time=[120., 145.])

sol_A03 = solve(prob_A03, alg=BOX(δ=0.04),
                clustering_method=LazyClustering(4),
                intersection_method=TemplateHullIntersection(boxdirs),
                intersect_source_invariant=false,
                intersect_source_invariant_method=TemplateHullIntersection(boxdirs),
                tspan = (0.0 .. 300.0))
property = SR02_specification(sol_A03)
push!(validation, Int(property))
SUITE[model][cases[7]] = @benchmarkable solve($prob_A03, alg=BOX(δ=0.04),
              clustering_method=LazyClustering(4),
              intersection_method=TemplateHullIntersection($boxdirs),
              intersect_source_invariant=false,
              intersect_source_invariant_method=TemplateHullIntersection($boxdirs),
              tspan = (0.0 .. 300.0))

# ----------------------------------------
#  A03 (discrete time)
# ----------------------------------------

sol_A03 = solve(prob_A03, alg=BOX(δ=0.1, approx_model=NoBloating()),
             clustering_method=LazyClustering(4),
             intersection_method=TemplateHullIntersection(boxdirs),
             intersect_source_invariant=false,
             intersect_source_invariant_method=TemplateHullIntersection(boxdirs),
             tspan = (0.0 .. 300.0))
property = SR02_specification(sol_A03)
push!(validation, Int(property))
SUITE[model][cases[8]] = @benchmarkable solve($prob_A03, alg=BOX(δ=0.1, approx_model=NoBloating()),
              clustering_method=LazyClustering(4),
              intersection_method=TemplateHullIntersection($boxdirs),
              intersect_source_invariant=false,
              intersect_source_invariant_method=TemplateHullIntersection($boxdirs),
              tspan = (0.0 .. 300.0))


# ----------------------------------------
#  A04 (dense time)
# ----------------------------------------

prob_A04 = spacecraft(abort_time=240.)

sol_A04 = solve(prob_A04, alg=BOX(δ=0.01),
              clustering_method=LazyClustering(40),
              intersection_method=TemplateHullIntersection(boxdirs),
              intersect_source_invariant=false,
              intersect_source_invariant_method=TemplateHullIntersection(boxdirs),
              tspan = (0.0 .. 300.0))
property = SR02_specification(sol_A04)
push!(validation, Int(property))
SUITE[model][cases[9]] = @benchmarkable solve($prob_A04, alg=BOX(δ=0.01),
            clustering_method=LazyClustering(40),
            intersection_method=TemplateHullIntersection($boxdirs),
            intersect_source_invariant=false,
            intersect_source_invariant_method=TemplateHullIntersection($boxdirs),
            tspan = (0.0 .. 300.0))

# ----------------------------------------
#  A04 (discrete time)
# ----------------------------------------

#=
sol_A04 = solve(prob_A04, alg=BOX(δ=0.1, approx_model=NoBloating()),
           clustering_method=LazyClustering(16),
           intersection_method=TemplateHullIntersection(boxdirs),
           intersect_source_invariant=false,
           intersect_source_invariant_method=TemplateHullIntersection(boxdirs),
           tspan = (0.0 .. 300.0))
property = SR02_specification(sol_A04)
push!(validation, Int(property))
SUITE[model][cases[10]] = @benchmarkable solve($prob_A04, alg=BOX(δ=0.1, approx_model=NoBloating()),
            clustering_method=LazyClustering(16),
            intersection_method=TemplateHullIntersection($boxdirs),
            intersect_source_invariant=false,
            intersect_source_invariant_method=TemplateHullIntersection($boxdirs),
            tspan = (0.0 .. 300.0))
=#

# ----------------------------------------
#  A05 (dense time)
# ----------------------------------------
#=
prob_A05 = spacecraft(abort_time=[235., 240.])

sol_A05 = solve(prob_A05, alg=BOX(δ=0.04),
              clustering_method=LazyClustering(13),
              intersection_method=TemplateHullIntersection(boxdirs),
              intersect_source_invariant=false,
              tspan = (0.0 .. 300.0))
property = SR02_specification(sol_A05)
push!(validation, Int(property))
SUITE[model][cases[11]] = @benchmarkable solve($prob_A05, alg=BOX(δ=0.04),
            clustering_method=LazyClustering(13),
            intersection_method=TemplateHullIntersection($boxdirs),
            intersect_source_invariant=false,
            tspan = (0.0 .. 300.0))
=#

# ----------------------------------------
#  A05 (discrete time)
# ----------------------------------------
#=
sol_A05 = solve(prob_A05, alg=BOX(δ=0.1, approx_model=NoBloating()),
           clustering_method=LazyClustering(13),
           intersection_method=TemplateHullIntersection(boxdirs),
           intersect_source_invariant=false,
           tspan = (0.0 .. 300.0))
property = SR02_specification(sol_A05)
push!(validation, Int(property))
SUITE[model][cases[12]] = @benchmarkable solve($prob_A05, alg=BOX(δ=0.1, approx_model=NoBloating()),
            clustering_method=LazyClustering(13),
            intersection_method=TemplateHullIntersection($boxdirs),
            intersect_source_invariant=false,
            tspan = (0.0 .. 300.0))
=#

# ----------------------------------------
#  A06 (dense time)
# ----------------------------------------
#=
prob_A06 = spacecraft(abort_time=[230., 240.])

sol_A06 = solve(prob_A06, alg=BOX(δ=0.04),
              clustering_method=LazyClustering(16),
              intersection_method=TemplateHullIntersection(boxdirs),
              intersect_source_invariant=false,
              tspan = (0.0 .. 300.0))
property = SR02_specification(sol_A06)
push!(validation, Int(property))
SUITE[model][cases[13]] = @benchmarkable solve($prob_A06, alg=BOX(δ=0.04),
            clustering_method=LazyClustering(16),
            intersection_method=TemplateHullIntersection($boxdirs),
            intersect_source_invariant=false,
            tspan = (0.0 .. 300.0))
=#

# ----------------------------------------
#  A06 (discrete time)
# ----------------------------------------
#=
sol_A06 = solve(prob_A06, alg=BOX(δ=0.1, approx_model=NoBloating()),
           clustering_method=LazyClustering(16),
           intersection_method=TemplateHullIntersection(boxdirs),
           intersect_source_invariant=false,
           tspan = (0.0 .. 300.0))
property = SR02_specification(sol_A06)
push!(validation, Int(property))
SUITE[model][cases[14]] = @benchmarkable solve($prob_A06, alg=BOX(δ=0.1, approx_model=NoBloating()),
            clustering_method=LazyClustering(16),
            intersection_method=TemplateHullIntersection($boxdirs),
            intersect_source_invariant=false,
            tspan = (0.0 .. 300.0))
=#

# ----------------------------------------
#  A07 (dense time)
# ----------------------------------------
#=
prob_A07 = spacecraft(abort_time=[50., 150.])

sol_A07 = solve(prob_A07, alg=BOX(δ=0.04),
              clustering_method=LazyClustering(50),
              intersection_method=TemplateHullIntersection(boxdirs),
              intersect_source_invariant=false,
              tspan = (0.0 .. 300.0))
property = SR02_specification(sol_A07)
push!(validation, Int(property))
SUITE[model][cases[15]] = @benchmarkable solve($prob_A07, alg=BOX(δ=0.04),
            clustering_method=LazyClustering(50),
            intersection_method=TemplateHullIntersection($boxdirs),
            intersect_source_invariant=false,
            tspan = (0.0 .. 300.0))
=#

# ----------------------------------------
#  A07 (discrete time)
# ----------------------------------------
#=
sol_A07 = solve(prob_A07, alg=BOX(δ=0.1, approx_model=NoBloating()),
           clustering_method=LazyClustering(50),
           intersection_method=TemplateHullIntersection(boxdirs),
           intersect_source_invariant=false,
           tspan = (0.0 .. 300.0))
property = SR02_specification(sol_A07)
push!(validation, Int(property))
SUITE[model][cases[16]] = @benchmarkable solve($prob_A07, alg=BOX(δ=0.1, approx_model=NoBloating()),
            clustering_method=LazyClustering(50),
            intersection_method=TemplateHullIntersection($boxdirs),
            intersect_source_invariant=false,
            tspan = (0.0 .. 300.0))
=#

# ----------------------------------------
#  A08 (dense time)
# ----------------------------------------
#=
prob_A08 = spacecraft(abort_time=[0., 240.])

sol_A08 = solve(prob_A08, alg=BOX(δ=0.04),
                clustering_method=LazyClustering(800),
                intersection_method=TemplateHullIntersection(boxdirs),
                intersect_source_invariant=false,
                tspan = (0.0 .. 300.0))
property = SR02_specification(sol_A08)
push!(validation, Int(property))
SUITE[model][cases[17]] = @benchmarkable solve($prob_A08, alg=BOX(δ=0.04),
            clustering_method=LazyClustering(800),
            intersection_method=TemplateHullIntersection($boxdirs),
            intersect_source_invariant=false,
            tspan = (0.0 .. 300.0))
=#

# ----------------------------------------
#  A08 (discrete time)
# ----------------------------------------
#=
sol_A08 = solve(prob_A08, alg=BOX(δ=0.1, approx_model=NoBloating()),
                clustering_method=LazyClustering(900),
                intersection_method=TemplateHullIntersection(boxdirs),
                intersect_source_invariant=false,
                tspan = (0.0 .. 300.0))
property = SR02_specification(sol_A08)
push!(validation, Int(property))
SUITE[model][cases[18]] = @benchmarkable solve($prob_A08, alg=BOX(δ=0.1, approx_model=NoBloating()),
            clustering_method=LazyClustering(900),
            intersection_method=TemplateHullIntersection($boxdirs),
            intersect_source_invariant=false,
            tspan = (0.0 .. 300.0))
=#

# ----------------------------------------
#  U01 (dense time)
# ----------------------------------------

prob_U01 = spacecraft(abort_time=260.)

sol_U01 = solve(prob_U01, alg=BOX(δ=0.04),
                clustering_method=LazyClustering(16),
                intersection_method=TemplateHullIntersection(boxdirs),
                intersect_source_invariant=false,
                tspan = (0.0 .. 300.0))
property = SR02_specification(sol_U01)
push!(validation, Int(property))
SUITE[model][cases[10]] = @benchmarkable solve($prob_U01, alg=BOX(δ=0.04),
            clustering_method=LazyClustering(16),
            intersection_method=TemplateHullIntersection($boxdirs),
            intersect_source_invariant=false,
            tspan = (0.0 .. 300.0))

# ----------------------------------------
#  U01 (discrete time)
# ----------------------------------------

#=
sol_U01 = solve(prob_U01, alg=BOX(δ=0.1, approx_model=NoBloating()),
                clustering_method=LazyClustering(16),
                intersection_method=TemplateHullIntersection(boxdirs),
                intersect_source_invariant=false,
                tspan = (0.0 .. 300.0))
property = SR02_specification(sol_U01)
push!(validation, Int(property))
SUITE[model][cases[12]] = @benchmarkable solve($prob_U01, alg=BOX(δ=0.1, approx_model=NoBloating()),
            clustering_method=LazyClustering(16),
            intersection_method=TemplateHullIntersection($boxdirs),
            intersect_source_invariant=false,
            tspan = (0.0 .. 300.0))
=#

# ----------------------------------------
#  U02 (dense time)
# ----------------------------------------

prob_U02 = spacecraft(abort_time=[0., 260.])

sol_U02 = solve(prob_U02, alg=BOX(δ=0.04),
                clustering_method=LazyClustering(16),
                intersection_method=TemplateHullIntersection(boxdirs),
                intersect_source_invariant=false,
                tspan = (0.0 .. 300.0))
property = SR02_specification(sol_U02)
push!(validation, Int(property))
SUITE[model][cases[11]] = @benchmarkable solve($prob_U02, alg=BOX(δ=0.04),
            clustering_method=LazyClustering(16),
            intersection_method=TemplateHullIntersection($boxdirs),
            intersect_source_invariant=false,
            tspan = (0.0 .. 300.0))

# ----------------------------------------
#  U02 (discrete time)
# ----------------------------------------

sol_U02 = solve(prob_U02, alg=BOX(δ=0.1, approx_model=NoBloating()),
                clustering_method=LazyClustering(16),
                intersection_method=TemplateHullIntersection(boxdirs),
                intersect_source_invariant=false,
                intersect_source_invariant_method=TemplateHullIntersection(boxdirs),
                tspan = (0.0 .. 300.0))
property = SR02_specification(sol_U02)
push!(validation, Int(property))
SUITE[model][cases[12]] = @benchmarkable solve($prob_U02, alg=BOX(δ=0.1, approx_model=NoBloating()),
            clustering_method=LazyClustering(16),
            intersection_method=TemplateHullIntersection($boxdirs),
            intersect_source_invariant=false,
            tspan = (0.0 .. 300.0))



sol_NA01 = nothing
sol_A01 = nothing
sol_A02 = nothing
sol_A03 = nothing
sol_A04 = nothing
#sol_A05 = nothing
#sol_A06 = nothing
#sol_A07 = nothing
#sol_A08 = nothing
sol_U01 = nothing
sol_U02 = nothing
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
    tm = median(results[model][c]).time * 1e-9
    runtimes[c] = round(tm, digits=4)
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

function plot_SRNA01()
    prob_NA01 = spacecraft(abort_time=-1.)
    sol = solve(prob_NA01, alg=BOX(δ=0.04),
                     clustering_method=LazyClustering(),
                     intersection_method=TemplateHullIntersection(boxdirs),
                     intersect_source_invariant=false,
                     tspan = (0.0 .. 300.0))

    idx_approaching = findall(x -> x == 1, location.(sol))
    idx_attempt = findall(x -> x == 2, location.(sol))
    idx_aborting = findall(x -> x == 3, location.(sol))

    fig = plot(legend=:bottomright, tickfont=font(30, "Times"), guidefontsize=45,
               xlab=L"x", ylab=L"y",
               xtick=[-800, -600, -400, -200.0, 0.0], ytick=[-400, -300, -200, -100, 0.],
               xlims=(-1000.0, 0.0), ylims=(-450.0, 0.0),
               bottom_margin=-6mm, left_margin=-3mm, right_margin=1mm, top_margin=3mm,
               size=(1000, 1000))

    for idx in idx_approaching
        plot!(fig, sol[idx], vars=(1, 2), lw=0.0, color=:lightgreen, alpha=1.)
    end
    for idx in idx_attempt
        plot!(fig, sol[idx], vars=(1, 2), lw=0.0, color=:red, alpha=1.)
    end
    for idx in idx_aborting
        plot!(fig, sol[idx], vars=(1, 2), lw=0.0, color=:cyan, alpha=1.)
    end
    fig
end

fig = plot_SRNA01()
savefig(fig, joinpath(@__DIR__, "ARCH-COMP22-JuliaReach-Spacecraft-SRNA01"))

# ----

function plot_SRA01()
    prob_A01 = spacecraft(abort_time=120.)
    @time sol = solve(prob_A01, alg=BOX(δ=0.04),
                      clustering_method=LazyClustering(3),
                      intersection_method=TemplateHullIntersection(boxdirs),
                      intersect_source_invariant=false,
                      tspan = (0.0 .. 300.0))

    idx_approaching = findall(x -> x == 1, location.(sol))
    idx_attempt = findall(x -> x == 2, location.(sol))
    idx_aborting = findall(x -> x == 3, location.(sol))

    fig = plot(legend=:bottomright, tickfont=font(30, "Times"), guidefontsize=45,
               xlab=L"x", ylab=L"y", lw=0.0,
               xtick=[-750, -500, -250, 0, 250.], ytick=[-400, -300, -200, -100, 0.],
               xlims=(-1000.0, 400.0), ylims=(-450.0, 0.0),
               bottom_margin=-6mm, left_margin=-3mm, right_margin=1mm, top_margin=3mm,
               size=(1000, 1000))

    for idx in idx_approaching
        plot!(fig, sol[idx], vars=(1, 2), lw=0.0, color=:lightgreen, alpha=1.)
    end
    for idx in idx_attempt
        plot!(fig, sol[idx], vars=(1, 2), lw=0.0, color=:red, alpha=1.)
    end
    for idx in idx_aborting
        plot!(fig, sol[idx], vars=(1, 2), lw=0.0, color=:cyan, alpha=1.)
    end
    fig
end

fig = plot_SRA01()
savefig(fig, joinpath(@__DIR__, "ARCH-COMP22-JuliaReach-Spacecraft-SRA01"))

GC.gc()
