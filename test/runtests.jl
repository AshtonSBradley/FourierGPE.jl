using FourierGPE, Test

@testset "Parseval test" begin include("test_parseval.jl") end
@testset "Ground state" begin include("test_groundstate.jl") end
@testset "Dynamics" begin include("test_dynamics.jl") end
@testset "Analysis" begin include("test_analysis.jl") end
