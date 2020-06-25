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

## total kinetic energy
k = LinRange(0.1,1.0,40) |> collect
Ek = kespectrum(k,ψi,X,K)

## autocorrelate ?
A = autocorrelate(ψi,X,K)

## ike spectrum
import FourierGPE:ikespectrum, kespectrum

function ikespectrum(k,ψ,X,K)
    x,y = X; kx,ky = K
 	dx,dy = diff(x)[1],diff(y)[1]
	DX,DK = dfftall(X,K)
    k2field = k2(K)
    psi = XField(ψ,X,K,k2field)
    vx,vy = velocity(psi)
    rho = abs2.(ψ)
    wx = @. sqrt(rho)*vx; wy = @. sqrt(rho)*vy
    Wi, Wc = helmholtz(wx,wy,psi)
    wix,wiy = Wi

	cwix = autocorrelate(wix,X,K)
	cwiy = autocorrelate(wiy,X,K)
    Ci = cwix .+ cwiy

    Nx = 2*length(x)
    Lx = x[end]-x[1] + dx
    xp = LinRange(-Lx,Lx,Nx+1)[1:Nx]
    yp = xp
    ρ = hypot.(xp,yp')

    Eki = zero(k)
    for i in eachindex(k)
        κ = k[i]
        Eki[i]  = 0.5*κ*sum(@. besselj0(κ*ρ)*Ci)*dx*dy |> real
    end
    return Eki
end

## call
Eki = ikespectrum(k,ψi,X,K)

## plot
plot(k,Eki)

##
ψi

## kespectrum fix
function kespectrum(k,ψ,X,K)
    x,y = X; kx,ky = K
    dx,dy = diff(x)[1],diff(y)[1]
    DX,DK = dfftall(X,K)
    k2field = k2(K)
    psi = XField(ψ,X,K,k2field)
    vx,vy = velocity(psi)
    rho = abs2.(ψ)
    wx = @. sqrt(rho)*vx; wy = @. sqrt(rho)*vy

    cwx = autocorrelate(wx,X,K)
    cwy = autocorrelate(wy,X,K)
    Ci = cwx .+ cwy

    Nx = 2*length(x)
    Lx = x[end]-x[1] + dx
    xp = LinRange(-Lx,Lx,Nx+1)[1:Nx]
    yp = xp
    ρ = hypot.(xp,yp')

    Eki = zero(k)
    for i in eachindex(k)
        κ = k[i]
        Eki[i]  = 0.5*κ*sum(@. besselj0(κ*ρ)*Ci)*dx*dy |> real
    end
    return Eki
end

## call
Ek = kespectrum(k,ψi,X,K)

## plot
plot(k,Ek)