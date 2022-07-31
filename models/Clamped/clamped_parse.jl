function parse_clamped_beam(model::String)
    aux = split(model, "-")
    N, Z = aux[2:3]
    N = parse(Int, N)
    path = "CB22$(Z)d_$(N).xml"
    H = readsxmodel(joinpath(@__DIR__, path), ST=CLCCS)
    mode = first(H.modes)
    return mode, N, Z
end
