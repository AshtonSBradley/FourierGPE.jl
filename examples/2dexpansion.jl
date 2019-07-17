using Pkg, Revise

#pkg"activate ."
using FourierGPE
gr(titlefontsize=12,size=(500,300),transpose=true,colorbar=false)

# ==== Units: ========================
# this example works in oscillator units

# ==== set simulation parameters ====
L = (40.0,40.0)
N = (256,256)
μ = 15.0

# ========= Initialize simulation ======
sim = Sim(L,N,par)
@pack! sim = μ
@unpack_Sim sim
# ===================================

# declare the potential function
import FourierGPE.V
V(x,y,t) = 0.5*(x^2 + y^2)

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

# turn off the potential
V(x,y,t) = 0.0


# ==== set simulation parameters ====
γ = 0.0
tf = 10.0
t = LinRange(ti,tf,Nt)
ϕi = ϕg

sime = Sim(L,N)
@pack! sime = tf,t,γ,ϕi
# initsim!(sim)

@unpack_Sim sime

# ====== Evolve in k space ==========
sole = runsim(sime)
# ===================================

ϕf = sole[end]
ψf = xspace(ϕf,sime)
showpsi(x,y,ψf)

#trim last few frames to show one orbit
# analytical result is within 10%
anim = @animate for i=1:Nt-6
    ψ = xspace(sole[i],sime)
    showpsi(x,y,ψ)
end

gif(anim,"./examples/expand.gif",fps=30)
