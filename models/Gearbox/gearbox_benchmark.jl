using BenchmarkTools, Plots, Plots.PlotMeasures, LaTeXStrings
import LazySets
using BenchmarkTools: minimum, median

SUITE = BenchmarkGroup()
model = "GEARBOX"
cases = ["GRBX01-MES01", "GRBX02-MES01", "GRBX01-MES01-discrete", "GRBX02-MES01-discrete"]

SUITE[model] = BenchmarkGroup()

include("gearbox.jl")
validation = []

LazySets.deactivate_assertions()

# Template directions
boxdirs = BoxDirections{Float64, Vector{Float64}}(6);
octdirs = CustomDirections([Vector(vi) for vi in OctDirections(6)]);

vx, vy, px, py, n, θ = 1, 2, 3, 4, 6, 0.628318530717959
vec1 = sparsevec([px, py], [tan(θ), 1.], n)
vec2 = sparsevec([px, py], [tan(θ), -1.], n)
vec1perp = sparsevec([px, py], [1, -tan(θ)], n)
vec2perp = sparsevec([px, py], [1., tan(θ)], n)
dirs_plane_34 = [vec1, -vec1, vec2, -vec2, vec1perp, -vec1perp, vec2perp, -vec2perp]

vec3 = sparsevec([vx, vy], [-sin(θ), -cos(θ)], n)
vec4 = sparsevec([vx, vy], [-sin(θ), cos(θ)], n)
vec3perp = sparsevec([vx, vy], [cos(θ), -sin(θ)], n)
vec4perp = sparsevec([vx, vy], [cos(θ), sin(θ)], n)

dirs_plane_12 = [vec3, -vec3, vec4, -vec4, vec3perp, -vec3perp, vec4perp, -vec4perp]
const extdirs = vcat(collect(boxdirs), Vector.(dirs_plane_34), Vector.(dirs_plane_12)) |> CustomDirections

# Initial states (non-homogeneous system)
X0_GRBX01 = Hyperrectangle(low=[0, 0, -0.0168, 0.0029, 0, 0], high=[0, 0, -0.0166, 0.0031, 0, 0])
X0_GRBX02 = Hyperrectangle(low=[0, 0, -0.01675, 0.00285, 0, 0], high=[0, 0, -0.01665, 0.00315, 0, 0]);

# Initial states (homogeneized system)
X0_GRBX01h = Hyperrectangle(low=[0, 0, -0.0168, 0.0029, 0, 1.], high=[0, 0, -0.0166, 0.0031, 0, 1.])
X0_GRBX02h = Hyperrectangle(low=[0, 0, -0.01675, 0.00285, 0, 1.], high=[0, 0, -0.01665, 0.00315, 0, 1.]);

thull = TemplateHullIntersection(extdirs);

# ----------------------------------------
#  GRBX01 (dense time)
# ----------------------------------------
prob_GRBX01 = gearbox_homog(X0=X0_GRBX01h)

# GRBX01-MES01
sol_GRBX01 = solve(prob_GRBX01, max_jumps=100,
                   clustering_method=LazyClustering(1, convex=false),
                   intersect_source_invariant=true,
                   intersection_source_invariant_method=thull,
                   intersection_method=thull,
                   tspan = 0 .. 0.21,
                   fixpoint_check=true,
                   alg=LGG09(δ=0.0005, template=extdirs, approx_model=Forward()))
property = req1(sol_GRBX01) && req2(sol_GRBX01)
push!(validation, Int(property))
SUITE[model][cases[1]] = @benchmarkable solve($prob_GRBX01, max_jumps=100,
                   clustering_method=LazyClustering(1, convex=false),
                   intersect_source_invariant=true,
                   intersection_source_invariant_method=$thull,
                   intersection_method=$thull,
                   tspan = 0 .. 0.21,
                   fixpoint_check=true,
                   alg=LGG09(δ=0.0005, template=$extdirs, approx_model=Forward()))

# ----------------------------------------
#  GRBX02 (dense time)
# ----------------------------------------
prob_GRBX02 = gearbox_homog(X0=X0_GRBX02h)

# GRBX02-MES01
sol_GRBX02 = solve(prob_GRBX02, max_jumps=100,
                   clustering_method=LazyClustering(1, convex=false),
                   intersect_source_invariant=true,
                   intersection_source_invariant_method=thull,
                   intersection_method=thull,
                   tspan = 0 .. 0.21,
                   fixpoint_check=true,
                   alg=LGG09(δ=0.0008, template=extdirs, approx_model=Forward()))
property = req1(sol_GRBX02) && req2(sol_GRBX02)
push!(validation, Int(property))
SUITE[model][cases[2]] = @benchmarkable solve($prob_GRBX02, max_jumps=100,
                   clustering_method=LazyClustering(1, convex=false),
                   intersect_source_invariant=true,
                   intersection_source_invariant_method=$thull,
                   intersection_method=$thull,
                   tspan = 0 .. 0.21,
                   fixpoint_check=true,
                   alg=LGG09(δ=0.0008, template=$extdirs, approx_model=Forward()))

# ----------------------------------------
#  GRBX01 (discrete time)
# ----------------------------------------
prob_GRBX01 = gearbox_homog(X0=X0_GRBX01h)

# GRBX01-MES01-discrete
sol_GRBX01d = solve(prob_GRBX01, max_jumps=100,
                   clustering_method=LazyClustering(1, convex=false),
                   intersect_source_invariant=true,
                   intersection_source_invariant_method=thull,
                   intersection_method=thull,
                   tspan = 0 .. 0.21,
                   fixpoint_check=true,
                   alg=LGG09(δ=0.0001, template=extdirs, cache=true, approx_model=NoBloating()))
property = req1(sol_GRBX01d) && req2(sol_GRBX01d)
push!(validation, Int(property))
SUITE[model][cases[3]] = @benchmarkable solve($prob_GRBX01, max_jumps=100,
                   clustering_method=LazyClustering(1, convex=false),
                   intersect_source_invariant=true,
                   intersection_source_invariant_method=$thull,
                   intersection_method=$thull,
                   tspan = 0 .. 0.21,
                   fixpoint_check=true,
                   alg=LGG09(δ=0.0001, template=extdirs, cache=true, approx_model=NoBloating()))

# ----------------------------------------
#  GRBX02 (discrete time)
# ----------------------------------------
prob_GRBX02 = gearbox_homog(X0=X0_GRBX02h)

# GRBX02-MES01-discrete
sol_GRBX02d = solve(prob_GRBX02, max_jumps=100,
                   clustering_method=LazyClustering(1, convex=false),
                   intersect_source_invariant=true,
                   intersection_source_invariant_method=thull,
                   intersection_method=thull,
                   tspan = 0 .. 0.21,
                   fixpoint_check=true,
                   alg=LGG09(δ=0.0001, template=extdirs, cache=true, approx_model=NoBloating()))
property = req1(sol_GRBX02d) && req2(sol_GRBX02d)
push!(validation, Int(property))
SUITE[model][cases[4]] = @benchmarkable solve($prob_GRBX02, max_jumps=100,
                   clustering_method=LazyClustering(1, convex=false),
                   intersect_source_invariant=true,
                   intersection_source_invariant_method=$thull,
                   intersection_method=$thull,
                   tspan = 0 .. 0.21,
                   fixpoint_check=true,
                   alg=LGG09(δ=0.0001, template=extdirs, cache=true, approx_model=NoBloating()))

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

fig = Plots.plot()
Plots.plot!(fig, sol_GRBX01, vars=(3, 4), ε=1e-5,
           color=:blue, alpha=0.5, lw=1.0, linecolor=:black,
           tickfont=font(30, "Times"), guidefontsize=45,
           xlab=L"x_3",
           ylab=L"x_4",
           xtick=([-0.016, -0.009, -0.002],
                  [L"-0.016", L"-0.009", L"-0.002"]),
           ytick=([-0.008, -0.004, 0.0, 0.004],
                  [L"-0.008", L"-0.004", L"0", L"0.004"]),
           xlims=(-0.017, -0.0015), ylims=(-0.008, 0.004),
           bottom_margin=-5mm, left_margin=-3mm, right_margin=10mm, top_margin=3mm,
           size=(1000, 1000))
savefig(fig, joinpath(@__DIR__, "ARCH-COMP22-JuliaReach-Gearbox_GRBX01.png"))

sol_GRBX01 = nothing
sol_GRBX02 = nothing
GC.gc()
