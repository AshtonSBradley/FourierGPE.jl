using Test, LinearAlgebra, BenchmarkTools, Parameters, Revise

N = 64
U0 = randn(N,N,N)

U1 = similar(U0)
U2 = similar(U0)

x = LinRange(-1,1,N)
y = x'

@time z = reshape(x,1,1,length(x))

V(x,y,z,t) = x^2 + y^2 + z^2

@time @. U0 += V(x,y,z,0.)

@time @. U0 = V(x,y,z,0.)

@time @. U1 += U0

@time @. U0 = U0 + U1

function Vadd!(V1,V0,V,t)
    @. V1 = V0 + V(x,y,z,t)
    return nothing
end

@time Vadd!(U1,U0,V,0.)

# function Valloc!(V1,V,t)
#     @. V1 = V(x,y,z,t)
#     return nothing
# end

# @btime Valloc!(V1,V,0.)

# macro

using FourierGPE
import FourierGPE.V

# static potential array
struct VStatic{N} <: UserParams
    V0::Array{Complex{Float64},N}
end

V0(x,y,z) = 0.25*(x^2 + y^2 + 4*z^2)
p0 = VStatic(V0.(x,y,z) .|> complex)


# dynamical potential function
@with_kw mutable struct Potential <: UserParams
    V::Expr = :( 0.5*(x^2 + y^2 + 4*z^2) )
end

macro potential(V,p)
    @eval $p = Potential($V)
end

macro staticpotential(V,p)
    @eval $p = VStatic($V.(x,y,z) .|> complex)
end

#TODO change sim to have V0, V, Vt fields
macro potentials(sim,p,q)
    @eval Sim($sim,V0 = $p,V = $q)
end

pot = Potential()
ex = pot.V

@potential :(0.25*(x^2 + y^2 + 4*z^2)) p1

@staticpotential V0 p2


#TODO V0 should just be a regular array
@test p2 == p0

p0 .≈ p2
