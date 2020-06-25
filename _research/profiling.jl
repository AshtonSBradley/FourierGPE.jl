# do some profiling
using Plots, LaTeXStrings
gr(colorbar=false,size=(600,150),legend=false,grid=false,xticks=false,yticks=false,axis=false)
msg =:mediumseagreen

using Revise, FourierGPE

# ==== Initialize simulation
L = (60.0,)
N = (512,)
sim = Sim(L,N)
@unpack_Sim sim

# ==== set simulation parameters
μ = 25.0
g = -0.01
γ = 0.0
Ns = 200
ξs = 2/abs(g)/Ns
us = 20
tf = π |> Float64
Nt = 150
ti = 0.0
t = LinRange(ti,tf,Nt)

# ==== soliton wavefunction
x = X[1]
ψs(x) = sqrt(Ns/2ξs)*sech(x/ξs)*exp(im*us*x)
ψi = ψs.(x)
ϕi = kspace(ψi,sim)

# ==== Set all fields
@pack_Sim! sim

# ====== Evolve in k space ==========
sol = runsim(sim)
# ===================================

import FourierGPE:Lgp!, nlin!, V
V(x,t) = 0.0

function nlin!(dϕ,ϕ,sim::Sim{1},t)
    @unpack g,X,V0 = sim; x = X[1]
    dϕ .= ϕ
    xspace!(dϕ,sim)
    @. dϕ *= V0 + V(x,[t]) + g*abs2(dϕ)
    kspace!(dϕ,sim)
    return nothing
end

function Lgp!(dϕ,ϕ,sim,t)
    @unpack γ,μ,espec = sim
    nlin!(dϕ,ϕ,sim,t)
    @. dϕ = -im*(1.0 - im*γ)*(dϕ + (espec - μ)*ϕ)
    return nothing
end

using BenchmarkTools

dϕi = copy(ϕi)
@btime Lgp!(dϕi,ϕi,sim,t)
# factor of 10 speed up in MKL in 1D
# if there isn't a speed up, there must be another bottleneck

@btime kspace!(dϕi,sim)
@btime xspace!(dϕi,sim)
@btime @. dϕi *= g*abs2(dϕ)


# ===============================
# 2d profiling: vortex precession example

using Revise, FourierGPE

# ==== Initialize simulation
L = (18.0,18.0)
N = (128,128)
sim = Sim(L,N)
@unpack_Sim sim

# ==== set simulation parameters
μ = 12.0

# Time dependent potential function (here trivial t dep)
import FourierGPE.V
V(x,y,t) = 0.5*(x^2 + y^2)
ψ0(x,y,μ,g) = sqrt(μ/g)*sqrt(max(1.0-V(x,y,0.0)/μ,0.0)+im*0.0)

# ==== make initial state
x,y = X
ψi = ψ0.(x,y',μ,g)
ϕi = kspace(ψi,sim)
@pack_Sim! sim

# ==== evolve
@time sol = runsim(sim)

# ==== ground state

ϕg = sol[end]
ψg = xspace(ϕg,sim)
showpsi(x,y,ψg)

# ==== set simulation parameters
γ = 0.0
t = LinRange(ti,tf,Nt)
ϕi = kspace(ψg,sim)
# reltol = 1e-7
# alg = DP5()ß

#---
using VortexDistributions
R(w) = sqrt(2*μ/w^2)
R(1)
rv = 3.
healinglength(x,y,μ,g) = 1/sqrt(g*abs2(ψ0(x,y,μ,g)))
ξ0 = healinglength.(0.,0.,μ,g)
ξ = healinglength(rv,0.,μ,g)
#---

pv = PointVortex(rv,0.,1)
vi = ScalarVortex(ξ,pv)
psi = Torus(copy(ψg),x,y)
vortex!(psi,vi)
showpsi(x,y,psi.ψ)
ψi .= psi.ψ
ϕi = kspace(ψi,sim)

# ===== compare with Fetter JLTP 2010
ξ = 1/sqrt(μ)
Rtf = R(1)
Ωm = 3*log(Rtf/ξ/sqrt(2))/2/Rtf^2
Ωv = Ωm/(1-rv^2/Rtf^2)

# or a precession period
Tv = 2*pi/Ωv

ti = 0.; tf = 0.1*Tv
t = LinRange(ti,tf,Nt)

@pack_Sim! sim
#---
# ==== evolve
@time solv = runsim(sim)


# FFTW 1.8s
# MKL 1.1s



import FourierGPE:Lgp!, nlin!, V
V(x,t) = 0.0

function nlin!(dϕ,ϕ,sim::Sim{1},t)
    @unpack g,X,V0 = sim; x = X[1]
    dϕ .= ϕ
    xspace!(dϕ,sim)
    @. dϕ *= V0 + V(x,[t]) + g*abs2(dϕ)
    kspace!(dϕ,sim)
    return nothing
end

function Lgp!(dϕ,ϕ,sim,t)
    @unpack γ,μ,espec = sim
    nlin!(dϕ,ϕ,sim,t)
    @. dϕ = -im*(1.0 - im*γ)*(dϕ + (espec - μ)*ϕ)
    return nothing
end

using BenchmarkTools

dϕi = copy(ϕi)
@btime Lgp!(dϕi,ϕi,sim,t)
# factor of ?? speed up in MKL in 2D
# if there isn't a speed up, there must be another bottleneck

@btime kspace!(dϕi,sim)
@btime xspace!(dϕi,sim)
@btime @. dϕi *= g*abs2(dϕ)
