##
using Test, SpecialFunctions, VortexDistributions
using FFTW, Plots, Pkg, Revise

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
rv = .8
healinglength(x,y,μ,g) = 1/sqrt(g*abs2(ψ0(x,y,μ,g)))
ξ0 = healinglength.(0.,0.,μ,g)
ξ = healinglength(rv,0.,μ,g)

pv = PointVortex(rv,0.,1)
nv = PointVortex(-rv,0,-1)
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

kmin = 0.1 #0.5*2*pi/R(1)
kmax = 2*pi/ξ0
Np = 200
kp = log10range(kmin,kmax,Np)

## find spec
Ek = kespectrum(kp,ψi,X,K)
plot(kp,abs.(Ek),scale=:log10,label="all KE",legend=:bottomleft)

Eki = ikespectrum(kp,ψi,X,K)
plot!(kp,abs.(Eki),scale=:log10,label="incompressible KE, nonperiodic")

Ekc = ckespectrum(kp,ψi,X,K)
plot!(kp,abs.(Ekc),scale=:log10,label="compressible KE")

# NOTE: spectra are not locally additive in k-space

## test totals (estimate since not linearly spaced)
dkp = diff(kp);push!(dkp,dkp[end])
Ekt = sum(Ek.*dkp)

## energy totals
K2 = k2(K)

function potential_energy(ϕ,sim,t)
    @unpack g,X = sim; x,y = X
    ψ = xspace(ϕ,sim)
    @. ψ *= V(x,y',t)
    return kspace(ψ,sim)
end

function interaction_energy(ϕ,sim)
    @unpack g,X = sim; x,y = X
    ψ = xspace(ϕ,sim)
    @. ψ *= 0.5*g*abs2(ψ)
    return kspace(ψ,sim)
end

function kinetic_energy(ϕ,sim)
    @unpack espec,K = sim; kx,ky = K
    dkx,dky = kx[2]-kx[1],ky[2]-ky[1]
    return sum(espec.*abs2.(ϕ))*dkx*dky |> real
end

dx = diff(x)[1]; dy = diff(y)[1]
dkx = diff(kx)[1]; dky = diff(ky)[1]
n0 = μ/g

function energies(ψ,sim)
    @unpack X,K = sim
    K2 = k2(K)
    psi = XField(ψ,X,K,K2)
    ek,ei,ec = energydecomp(psi)
    Ekh = sum(ek)*dx*dy |> real
    Ei = sum(ei)*dx*dy |> real
    Ec = sum(ec)*dx*dy |> real
    Ek = kinetic_energy(kspace(ψ,sim),sim)

    Ev = sum(@. V(x,y',t)*abs2(ψ))*dx*dy |> real
    Eint = sum(@. g/2*(abs2(ψ) - n0)^2)*dx*dy |> real
    Eqp = Ekall - Ekhy

    Et = Ekall + Ev + Eint
    Natoms = sum(abs2.(ψ))*dx*dy
    # end
    return Ei,Ec,Ekhy,Eqp,Ev,Eint,Et,Ekall,Natoms
end

Ei,Ec,Ekhy,Eqp,Ev,Eint,Et,Ekall,Natoms = energies(psi.ψ,sim)

##
ϕ = kspace(psi.ψ,sim)
ke = kinetic_energy(kspace(psi.ψ,sim),sim)

## check ordering with plot
using Plots
heatmap(log10.(abs2.(fftshift(ϕ))))
##
heatmap(real(abs.(fftshift(sim.espec))))
# OK, so everything in kspace is in fft order

##
