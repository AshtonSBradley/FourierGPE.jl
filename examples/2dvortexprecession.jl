using Revise, FourierGPE, ColorSchemes
gr(titlefontsize=12,size=(500,300),transpose=true,colorbar=false)
turbo = cgrad(ColorSchemes.turbo.colors)
import FourierGPE.showpsi

function showpsi(x,y,ψ)
    p1 = heatmap(x,y,abs2.(ψ),aspectratio=1,c=turbo)
    xlims!(x[1],x[end]);ylims!(y[1],y[end])
    xlabel!(L"x");ylabel!(L"y")
    title!(L"|\psi|^2")
    p2 = heatmap(x,y,angle.(ψ),aspectratio=1,c=turbo)
    xlims!(x[1],x[end]);ylims!(y[1],y[end])
    xlabel!(L"x");ylabel!(L"y")
    title!(L"\textrm{phase} (\psi)")
    p = plot(p1,p2,size=(600,300))
    return p
end

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

vcore = Exact(ξ)
pv = PointVortex(rv,0.,1)
vi = ScalarVortex(vcore,pv)
psi = Torus(copy(ψg),x,y)
vortex!(psi,vi)
showpsi(x,y,psi.ψ)
ψi .= psi.ψ
ϕi = kspace(ψi,sim)
showpsi(x,y,psi.ψ)

# diploe
# pv = PointVortex(rv,0.,1)
# nv = PointVortex(-rv,0.,-1)
# dipole = [pv;nv]
# di = ScalarVortex.([vcore],dipole) #TODO missing method!
# psi = Torus(copy(ψg),x,y)
# vortex!(psi,di)

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
solv = runsim(sim)

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
