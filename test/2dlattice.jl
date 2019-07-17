using Plots, LaTeXStrings, Pkg, Revise
gr(titlefontsize=12,size=(500,300),transpose=true,colorbar=false)

#pkg"activate ."
using FourierGPE

# ==== Units: ========================
# this example works in oscillator units


# ==== set simulation parameters ====
L = (20.0,20.0)
N = (128,128)
μ = 25.0

# ========= Initialize simulation ======
sim = Sim(L,N)
@pack! sim = μ
@unpack_Sim sim


# declare the potential function
import FourierGPE.V
V(x,y,t) = 0.5*(x^2 + y^2)

# useful TF state
using VortexDistributions
ψ0(x,y,μ,g) = sqrt(μ/g)*sqrt(max(1.0-V(x,y,0.0)/μ,0.0)+im*0.0)

#make initial state
x,y = X
ψi = ψ0.(x,y',μ,g)
makevortex!(ψi,[0. 0. 1],x,y)
showpsi(x,y,ψi)
ϕi = kspace(ψi,sim)
@pack! sim = ϕi
sim
@unpack T = sim

T = Transforms(T,psi = ArrayPartition(ϕi,copy(ϕi),copy(ϕi)))

sim = Sim(sim,T=T)
ϕj = copy(ϕi)

function Lz(ϕ,sim)
    @unpack T,X,K = sim
    @unpack psi = T
    xky,ykx,L = psi.x[1],psi.x[2],psi.x[3]
    kx,ky = K
    x,y = X
    xky .= ϕ.*ky'
    ykx .= kx.*ϕ
    T.Tkx!*xky
    T.Tkx!*ykx
    @. L = x*xky - ykx*y'
    T.Txk!*L
    return L
end

@time Ltest = Lz(ϕj,sim)
Ltest == ϕj
ψjL = xspace(Ltest,sim)
sum(conj(ψi).*ψjL)/sum(abs2.(ψi))

# Ok, we have allocation free evaluation!!!

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
