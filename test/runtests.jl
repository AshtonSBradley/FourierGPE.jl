using FourierGPE, Test

#runtests
@testset "Transform tests" begin include("transformtests.jl") end
# @testset "Ground state" begin include("groundstatetests.jl") end
# @testset "Evolution tests" begin include("evolutiontests.jl") end
