using FourierGPE, Test

#runtests
@testset "Parseval tests" begin include("parsevaltests.jl") end
# @testset "Mixed Parseval tests" begin include("mixedparsevaltests.jl") end
@testset "Ground state tests" begin include("groundstatetests.jl") end
