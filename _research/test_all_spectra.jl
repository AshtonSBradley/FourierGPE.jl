## load
using Test, SpecialFunctions, VortexDistributions
using FourierGPE, Plots

## Initialize simulation
L = (18.0,18.0)
N = (256,256)
sim = Sim(L,N)
@unpack_Sim sim

## make initial state
x,y = X
ψi = 100*one.(x.*y') |> complex
ϕi = kspace(ψi,sim)
@pack_Sim! sim
dx,dy = diff(x)[1],diff(y)[1]

## plot
heatmap(x,y,abs2.(ψi))


## call
k = LinRange(0.01,10.0,100) |> collect
Eki = ikespectrum(k,ψi,X,K)

## plot
plot(k,Eki)



## call

Ek = kespectrum(k,ψi,X,K)

## plot
plot(k,Ek)

## what does autocorrelate look like?
x,y = X; kx,ky = K
dx,dy = diff(x)[1],diff(y)[1]
psiac = autocorrelate(ψi,X,K)
psiac2 = autocorrelate(ψi,X,K,periodic=true)

## plot
heatmap(abs2.(psiac2))

##
psiac2