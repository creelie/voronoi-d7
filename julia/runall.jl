for f in ("E7Cell.jl", "RootLattices.jl", "SolidAngle.jl")
    println("== ", f); include(joinpath(@__DIR__, f))
end
