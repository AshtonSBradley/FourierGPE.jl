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
    for (j,kj) in enumerate(k)
        Eki[j]  = 0.5*kj*sum(@. besselj0(kj*ρ)*Ci)*dx*dy |> real
    end
    return Eki
end

## call
Eki = ikespectrum(k,ψi,X,K)

## plot
plot(k,Eki)

## kespectrum fix
function kespectrum(k,ψ,X,K)
    x,y = X; kx,ky = K
    dx,dy = diff(x)[1],diff(y)[1]
    psiac = autocorrelate(ψ,X,K,periodic=true)

    # make ρ
    Nx = 2*length(x)
    Lx = x[end]-x[1] + dx
    xp = LinRange(-Lx,Lx,Nx+1)[1:Nx]
    yp = xp
    ρ = hypot.(xp,yp')

    # do spectra
    Ek = zero(k)
    for (i,ki) in enumerate(k)
        Ek[i]  = 0.5*ki^3*sum(@. besselj0(ki*ρ)*psiac)*dx*dy |> real
    end
    return Ek
end

## call
k = LinRange(0.01,10.0,100) |> collect
Ek = kespectrum(k,ψi,X,K)

## plot
plot(k,Ek)

## make autocorrelate periodic
using FFTW
import FourierGPE:autocorrelate

@doc raw"""
	autocorrelate(ψ,X,K,periodic=false)

Return the autocorrelation integral of a complex field ``\psi``, ``A``, given by

```
A(\rho)=\int d^2r\;\psi^*(r-\rho)\psi(r)
```

defined on a cartesian grid on a cartesian grid using FFTW.

`X` and `K` are tuples of vectors `x`,`y`,`kx`, `ky`.

This method is useful for evaluating spectra from cartesian data.
"""
function autocorrelate(ψ,X,K;periodic=false)
    DX,DK = dfftall(X,K)
    if periodic == false
        ϕ = zeropad(ψ)
    else
        ϕ = ψ
    end
	χ = fft(ϕ)*prod(DX)
	return ifft(abs2.(χ))*prod(DK) |> fftshift
end

## what does autocorrelate look like?
x,y = X; kx,ky = K
dx,dy = diff(x)[1],diff(y)[1]
psiac = autocorrelate(ψi,X,K)
psiac2 = autocorrelate(ψi,X,K,periodic=true)

## plot
heatmap(abs2.(psiac2))

##
psiac2