using FFTW, Parameters

abstract type Field end

struct XField{D} <: Field
    psiX::Array{Complex{Float64},D}
    X::NTuple{D}
    K::NTuple{D}
    K2::Array{Float64,D}
end

struct KField{D} <: Field
    psiK::Array{Complex{Float64},D}
    X::NTuple{D}
    K::NTuple{D}
    K2::Array{Float64,D}
end

xvec(L,N) = LinRange(-L/2,L/2,N+1)[2:end] |> collect

function kvec(λ,N)
    @assert iseven(N)
    nkx = 0:Int(N/2)
    kx = [nkx[1:end-1];-reverse(nkx[2:end])]*2*π/λ
return kx
end

function xvecs(L,N)
    X = []
    for (λ,ν) in zip(L,N)
        x = xvec(λ,ν)
        push!(X,x)
    end
    return X |> Tuple
end

function kvecs(L,N)
    K = []
    for (λ,ν) in zip(L,N)
        k = kvec(λ,ν)
        push!(K,k)
    end
    return K |> Tuple
end

function k2(K)
    kind = Iterators.product(K...)
    return map(x-> sum(abs2.(x)),kind)
end

include("../src/analysis.jl")

N = 100

X = xvecs((1,1),(N,N))
K = kvecs((1,1),(N,N))
K2 = k2(K)

ktest = 2*pi
@. psi = exp(im*ktest*X[1]*one.(X[2]'))
# Note: ktest = n*2*pi that are evaluated exactly are the precise values of representation.
# Hence any derivative of a field constructed from a superpoisition of these k's will
# also be exact.

#2d test
xpsi = XField(psi,X,K,K2)
kpsi = KField(fft(psi),X,K,K2)

vx,vy = velocity(xpsi)

@time vx,vy = velocity(xpsi)

#3d test
X = xvecs((1,2,3),(N,N,N))
K = kvecs((1,2,3),(N,N,N))
K2 = k2(K)
psi = randn(N,N,N) + im*randn(N,N,N)

xpsi = XField(psi,X,K,K2)
kpsi = KField(fft(psi),X,K,K2)

@time vx,vy,vz = velocity(xpsi)
