using LinearAlgebra, BenchmarkTools, Revise

N = 512
V0 = randn(N,N,N)

V1 = similar(V0)
V2 = similar(V0)

x = LinRange(-1,1,N)
y = x'

z = reshape(x,1,1,length(x))

Vt(x,y,z) = x^2 + y^2 + z^2

@btime @. V0 += Vt(x,y,z)

@btime @. V0 = Vt(x,y,z)

@btime @. V1 += V0

@btime @. V0 = V0 + V1

function addpotentials(V0,V1)
    @inbounds for i in eachindex(V0)
        V0[i] = V0[i] + V1[i]
    end
    return V0
end

@btime V2 = addpotentials(V0,V1)

function makepotential(V0,Vt,x,y,z)
    @. V0 += Vt(x,y,z)
    return V0
end

@btime V2 = makepotential(V0,Vt,x,y,z)


# macro

using FourierGPE
import FourierGPE.V

@with_kw mutable struct Potential <: UserParams @deftype Float64
    V::Expr = :( V(x,y,z,t) = 0.5*(x^2 + y^2 + 4*z^2) )
end

pot = Potential()
# Symbol(par.V)
# import FourierGPE.V
# eval(par.V)

macro potential(V)
    @eval p = Potential($V)
    @eval $p.V
    return p
end

@potential :(V(x,y,z,t) = 0.5*(x^2 + y^2 + 4*z^2))


macro capture(fdef)
    @eval $fdef
end

@capture :(f(x)=x^2)

f(1)
