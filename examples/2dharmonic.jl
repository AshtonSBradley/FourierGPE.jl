using Plots, LaTeXStrings, Pkg, Revise
gr(titlefontsize=12,size=(500,300),transpose=true,colorbar=false)

#pkg"activate ."
using FourierGPE

# ==== Units: ========================
# this example works in oscillator units
# convenient plot
function showpsi(x,y,ψ)
    p1 = heatmap(x,y,abs2.(ψ),aspectratio=1)
    xlabel!(L"x/a_x");ylabel!(L"y/a_y")
    title!(L"|\psi|^2")
    p2 = heatmap(x,y,angle.(ψ),aspectratio=1)
    xlabel!(L"x/a_x");ylabel!(L"y/a_y")
    title!(L"\textrm{phase}(\psi)")
    p = plot(p1,p2,size=(600,300))
    return p
end

# ==== define user parameters =======
@with_kw mutable struct Params <: UserParams @deftype Float64
    # user parameters:
    κ = 0.1
end
par = Params()

# ==== set simulation parameters ====
L = (20.0,20.0)
N = (128,128)
μ = 25.0
# X,K,dX,dK,DX,DK,T = maketransforms(L,N)
# espec = 0.5*k2(L...,N...)

# ===# ====== Initialize simulation ======
sim = Sim(L,N,par)
@pack! sim = μ
@unpack_Sim sim
# ====================================== Initialize simulation ======
# sim = Sim(L,N,par)
# @pack! sim = T,X,K,espec
# initsim!(sim)
# @unpack_Sim sim
# ===================================

# declare the potential function
import FourierGPE.V
V(x,y,t)::Float64 = 0.5*(x^2 + y^2)

# useful TF state
ψ0(x,y,μ,g) = sqrt(μ/g)*sqrt(max(1.0-V(x,y,0.0)/μ,0.0)+im*0.0)

x,y = X
#make initial state
ψi = ψ0.(x,y',μ,g)
ϕi = kspace(ψi,sim)
@pack! sim = ϕi
sim

# ====== Evolve in k space ==========
sol = runsim(sim)
# ===================================

# pull out the ground state:
ϕg = sol[end]
ψg = xspace(ϕg,sim)
showpsi(x,y,ψg)

# Add a vortex off-axis
using VortexDistributions

# initial state: imprint a vortex inside Thomas-Fermi radius
healing(x,y,μ,g) = 1/sqrt(g*abs2(ψ0(x,y,μ,g)))
Rtf = sqrt(2*μ)
rv = 0.5*Rtf
xv,yv,cv = rv, 0.0, 1
initialvortex = [xv yv cv]
ξv = healing(xv,yv,μ,g)
ψv = copy(ψg)
makeallvortices!(ψv,initialvortex,x,y,ξv)
showpsi(x,y,ψv)

# In TF regim precession frequency is given analytically by:
# (see Fetter JLTP 2010)
ξ = 1/sqrt(μ)
Ωm = 3*log(Rtf/ξ/sqrt(2))/2/Rtf^2
Ωv = Ωm/(1-rv^2/Rtf^2)

#or a period of
Tv = 2*π/Ωv

# ==== set simulation parameters ====
γ = 0.0
tf = Tv
t = LinRange(ti,tf,Nt)
ϕi = kspace(ψv,sim)
reltol = 1e-7
alg = DP5()
@pack! sim = tf,t,γ,ϕi,reltol,alg
# initsim!(sim)

@unpack_Sim sim

# ====== Evolve in k space ==========
solv = runsim(sim)
# ===================================

ϕf = solv[100]
ψf = xspace(ϕf,sim)
showpsi(x,y,ψf)

#trim last few frames to show one orbit
# analytical result is within 10%
anim = @animate for i=1:Nt-6
    ψ = xspace(solv[i],sim)
    showpsi(x,y,ψ)
end

gif(anim,"./examples/vortex.gif",fps=30)
