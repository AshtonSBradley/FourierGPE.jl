using FourierGPE, Test

@testset "Parseval test" begin include("parseval_tests.jl") end
@testset "Mixed Parseval" begin include("mixedparseval_tests.jl") end
@testset "Ground state" begin include("groundstate_tests.jl") end
@testset "Dynamics" begin include("dynamics_tests.jl") end
@testset "Analysis" begin include("analysis_tests.jl") end
