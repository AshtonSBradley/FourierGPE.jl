using FourierGPE, Test

@testset "Parseval test" begin include("parsevaltests.jl") end
# @testset "Mixed Parseval" begin include("mixedparsevaltests.jl") end
@testset "Ground state" begin include("groundstatetests.jl") end
@testset "Dynamics" begin include("dynamicstests.jl") end
