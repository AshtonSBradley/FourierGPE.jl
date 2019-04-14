#transform tests
using Revise

using FourierGPE

function makearrays(L,N)
    X = xvecs(L,N)
    K = kvecs(L,N)
    dX = Float64[]; dK = Float64[]
    for j ∈ eachindex(X)
        x=X[j];k=K[j]
        dx = x[2]-x[1];dk = k[2]-k[1]
        push!(dX,dx)
        push!(dK,dk)
    end
    dX = dX |> Tuple
    dK = dK |> Tuple
    return X,K,dX,dK
end

L = (11.,140.)
N = (20,40)
X,K,dX,dK = makearrays(L,N)
DX,DK = dfftall(X,K)

#measures for unitary FFT
dμx = prod(DX); dμk = prod(DK)

# plan transforms
FFTW.set_num_threads(Sys.CPU_THREADS)
flags = FFTW.MEASURE
ψtest = ones(N...) |> complex

Tij       = [:Txk :Txk! :Tkx :Tkx!]
dμi       = [:dμx :dμx :dμk :dμk]
trans     = [:fft :fft! :ifft :ifft!]
prefix    = :plan_

for (Tij,dμi,trans) ∈ zip(Tij,dμi,trans)
    plantrans = Symbol(prefix,trans)
    @eval $Tij = $dμi*$plantrans(ψtest,flags=flags)
end

T = Transforms(Txk,Txk!,Tkx,Tkx!)

# Txk = dμx*plan_fft(ψtest,flags=flags)
# # ψtest = ones(Ns...) |> complex
# Txk! = dμx*plan_fft!(ψtest,flags=flags)
# # ψtest = ones(Ns...) |> complex
# Tkx  = dμk*plan_ifft(ψtest,flags=flags)
# # ψtest = ones(Ns...) |> complex
# Tkx!  = dμk*plan_ifft!(ψtest,flags=flags)
# # ψtest = ones(Ns...) |> complex
