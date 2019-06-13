using FFTW

abstract type Field end

struct XField{D} <: Field
    XPsi::Array{Complex{Float64},D}
    XGrids::NTuple{D}
end

struct KField{D} <: Field
    KPsi::Array{Complex{Float64},D}
    KGrids::NTuple{D}
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

N = 100
psi = randn(N,N) + im*randn(N,N)

X = xvecs((1,2),(N,N))
K = kvecs((1,2),(N,N))
K2 = k2(K)

typeof(K)

xpsi = XField(psi,X)
kpsi = KField(fft(psi),K,K2)

sizeof(xpsi.XPsi)

sizeof(psi)
