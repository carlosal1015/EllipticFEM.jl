using EllipticFEM
using Base.Test

@test include("./Elliptic/runElliptic.jl") == nothing
