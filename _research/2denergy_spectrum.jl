##
using Test, SpecialFunctions, VortexDistributions
using FFTW, Plots, Pkg, LaTeXStrings

##
Pkg.activate(".")
using FourierGPE

## Initialize simulation
# harmonic oscillator units
L = (18.0,18.0)
N = (256,256)
sim = Sim(L,N)
@unpack_Sim sim

# set simulation parameters
μ = 20.0

## Time dependent potential function (here trivial t dep)
import FourierGPE.V
V(x,y,t) = 0.5*(x^2 + y^2)
ψ0(x,y,μ,g) = sqrt(μ/g)*sqrt(max(1.0-V(x,y,0.0)/μ,0.0)+im*0.0)

## make initial state
x,y = X
ψi = ψ0.(x,y',μ,g)
ϕi = kspace(ψi,sim)
@pack_Sim! sim
dx,dy = diff(x)[1],diff(y)[1]

## evolve
@time sol = runsim(sim)

## pull out ground state
ϕg = sol[end]
ψg = xspace(ϕg,sim)
showpsi(x,y,ψg)

## make dipole
R(w) = sqrt(2*μ/w^2)
R(1)
d = 2
healinglength(x,y,μ,g) = 1/sqrt(g*abs2(ψ0(x,y,μ,g)))
ξ0 = healinglength.(0.,0.,μ,g)
ξ = healinglength(d/2,0.,μ,g)

pv = PointVortex(d/2,0.,1)
nv = PointVortex(-d/2,0,-1)
v1 = ScalarVortex(ξ,pv)
v2 = ScalarVortex(ξ,nv)

psi = Torus(copy(ψg),x,y)
vortex!(psi,[v1;v2])
showpsi(x,y,psi.ψ)

## kespectrum
ψi .= psi.ψ
ϕi = kspace(ψi,sim) |> fftshift
kx,ky = K .|> fftshift
x,y = X

kR = 2*pi/R(1)
kxi = 2*pi/ξ0

kmin = 0.15kR
kmax = kxi
Np = 200
kp = log10range(kmin,kmax,Np)

## find spec
Ek = kespectrum(kp,ψi,X,K)
plot(kp*ξ0,Ek,scale=:log10,label=L"E_{kin}(k)",legend=:bottomleft,grid=false,foreground_color_legend = nothing)

Eki = ikespectrum(kp,ψi,X,K)
plot!(kp*ξ0,Eki,scale=:log10,label=L"\epsilon_k^i(k)")

Ekc = ckespectrum(kp,ψi,X,K)
plot!(kp*ξ0,Ekc,scale=:log10,label=L"\epsilon_k^c(k)")

xlabel!(L"k\xi")

# k values of interest
vline!([(2*pi/R(1))*ξ0],label=L"2\pi/R")
vline!([(2*pi/d)*ξ0],label=L"2\pi/d")
vline!([1],ls=:dash,label=L"1/\xi")
# vline!([(2*pi/dx)*ξ0],label=L"2\pi/dx") # check grid cutoff


# Eqp = qpespectrum(kp,ψi,X,K)
# plot!(kp,Eqp,scale=:log10,label="quantum pressure")
# NOTE: spectra are not locally additive in k-space
# Testing if decomp v versus u: v decomp generates a lot of spurious energy at k->0? And compressible energy would seem to be enormous for an incompressible state. No pathology, but incorrect measure. 

## compute total kinetic energy spectrally from momentum space w.f. ϕ
function kinetic_energy(ϕ,sim)
    @unpack espec,K = sim; kx,ky = K
    dkx,dky = kx[2]-kx[1],ky[2]-ky[1]
    return sum(espec.*abs2.(ϕ))*dkx*dky |> real
end

## test totals against spectral evauluation
dkp = diff(kp);push!(dkp,dkp[end])

# integrate each
Ekit = sum(Eki.*dkp)
Ekct = sum(Ekc.*dkp)
Eqpt = sum(Eqp.*dkp)
#sum
Ekit + Ekct + Eqpt

#integrate total spectrum
Ekt = sum(Ek.*dkp)

#integrate total in k space
ke = kinetic_energy(kspace(psi.ψ,sim),sim)




## NOTE stop here

# some convenience methods
# dx = diff(x)[1]; dy = diff(y)[1]
# dkx = diff(kx)[1]; dky = diff(ky)[1]
# n0 = μ/g

function energies(ψ,sim)
    @unpack X,K = sim
    x,y = X
    K2 = k2(K)
    psi = XField(ψ,X,K,K2)
    ek,ei,ec = energydecomp(psi)
    Ekh = sum(ek)*dx*dy |> real
    Ei = sum(ei)*dx*dy |> real
    Ec = sum(ec)*dx*dy |> real
    Ek = kinetic_energy(kspace(ψ,sim),sim)
    # Ev = sum(@. V(x,y',t)*abs2(ψ))*dx*dy |> real
    Eint = sum(@. g/2*(abs2(ψ) - n0)^2)*dx*dy |> real
    Eqp = Ek - Ekh

    Et = Ek + Eint # + Ev
    Natoms = sum(abs2.(ψ))*dx*dy
    # end
    return Ei,Ec,Ekh,Eqp,Eint,Et,Ek,Natoms
end

## evaluate
Ei,Ec,Ekh,Eqp,Eint,Et,Ek,Natoms = energies(psi.ψ,sim)
