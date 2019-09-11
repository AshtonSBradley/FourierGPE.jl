using Revise, FourierGPE

# ==== Initialize simulation
L = (25.0,25.0)
N = (256,256)
sim = Sim(L,N)
@unpack_Sim sim

# ==== set simulation parameters
μ = 12.0

# Time dependent potential function (here trivial t dep)
import FourierGPE.V
V(x,y,t) = 0.5*(x^2 + y^2)
ψ0(x,y,μ,g) = sqrt(μ/g)*sqrt(max(1.0-V(x,y,0.0)/μ,0.0)+im*0.0)

# Or, as a static Potential
# x,y = X
# Vs(x,y) = 0.5*(x^2 + y^2)
# ψ0(x,y,μ,g) = sqrt(μ/g)*sqrt(max(1.0-Vs(x,y)/μ,0.0)+im*0.0)
# sim = Sim(sim,V0=Vs.(x,y'))

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
# alg = DP5()

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
showpsi(x,y,psi.ψ)

# import VortexDistributions:ScalarVortex,vortex!
#
# ScalarVortex(ξ::Float64,pv) = ScalarVortex.([Exact(ξ)],pv)
# ScalarVortex(ξ::Array{Float64,1},pv) = @. ScalarVortex(Exact(ξ),pv)
#
# ξa = [ξ;2*ξ]
#
# Exact.(ξa)
# # dipole test
# pv = PointVortex(rv,0.,1)
# nv = PointVortex(-rv,0.,-1)
# dipole = [pv;nv]
#
# ScalarVortex(ξ,dipole)
# ScalarVortex(ξa,dipole)

# di = ScalarVortex(ξ,dipole) #TODO missing method!
# di = ScalarVortex([ξ1;ξ2],dipole)
# vortex!(psi <: Array{Complex{Float64},2},nx3 array of coordinates) return an array

# di = ScalarVortex.([vcore(ξ)],dipole) #TODO missing method!
# di = ScalarVortex.([vcore(ξ1);vcore(ξ2)],dipole)

# psi = Torus(copy(ψg),x,y)
# vortex!(psi,di)
showpsi(x,y,psi.ψ)
#---
# ===== compare with Fetter JLTP 2010
ξ = 1/sqrt(μ)
Rtf = R(1)
Ωm = 3*log(Rtf/ξ/sqrt(2))/2/Rtf^2
Ωv = Ωm/(1-rv^2/Rtf^2)

# or a precession period
Tv = 2*pi/Ωv

ti = 0.; tf = Tv
t = LinRange(ti,tf,Nt)

@pack_Sim! sim
#---
# ==== evolve
@time solv = runsim(sim)

# ==== analyse
ϕf = solv[end]
ψf = xspace(ϕf,sim)
showpsi(x,y,ψf)

anim = @animate for i in eachindex(t)
    ψ = xspace(solv[i],sim)
    showpsi(x,y,ψ)
end

saveto = joinpath(@__DIR__,"vortexprecession.gif")
gif(anim,saveto,fps=25)
