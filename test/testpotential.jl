using LinearAlgebra, BenchmarkTools

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
