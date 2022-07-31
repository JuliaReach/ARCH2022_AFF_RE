# deactivate plot GUI, which is not available in Docker
ENV["GKSwstype"] = "100"

# instantiate project
import Pkg
Pkg.activate(@__DIR__)
Pkg.instantiate()

const TEST_LONG = true   # if true, the longer test suite is run; may take > 1 hour
                         # and requires at least 16gb RAM
const TARGET_FOLDER = "result"
const RESULTS_FILE = "results.csv"

function main()
    if !isdir(TARGET_FOLDER)
        mkdir(TARGET_FOLDER)
    end
    global io = open(joinpath(TARGET_FOLDER, RESULTS_FILE), "w")

    println("Running AFF benchmarks...")

    # Heat 3D benchmark
    println("###\nRunning Heat 3D benchmark\n###")
    include("models/Heat3D/heat3d_benchmark.jl")

    # ISS benchmark
    println("###\nRunning ISS benchmark\n###")
    include("models/ISS/ISS_benchmark.jl")

    # Gearbox benchmark
    println("###\nRunning Gearbox benchmark\n###")
    include("models/Gearbox/gearbox_benchmark.jl")

    # Spacecraft benchmark
    if TEST_LONG
        println("###\nRunning spacecraft benchmark\n###")
        include("models/Spacecraft/Spacecraft_benchmark.jl")
    end

    # Building benchmark
    println("###\nRunning building benchmark\n###")
    include("models/Building/building_benchmark.jl")

    # Platoon benchmark
    println("###\nRunning platoon benchmark\n###")
    include("models/Platoon/Platoon_benchmark.jl")

    # Electromechanic break benchmark
    if TEST_LONG
        println("###\nRunning electromechanic-break benchmark\n###")
        include("models/EMBrake/embrake_benchmark.jl")
    end

    # Powertrain benchmark
    println("###\nRunning Powertrain benchmark\n###")
    include("models/Powertrain/Powertrain_benchmark.jl")

    # Clamped benchmark
    println("###\nRunning Clamped benchmark\n###")
    include("models/Clamped/clamped_benchmark.jl")

    print(io, "\n")
    println("Finished running benchmarks.")
    close(io)
    nothing
end

main()
