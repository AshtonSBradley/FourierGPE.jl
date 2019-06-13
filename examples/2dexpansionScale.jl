using Pkg, Revise

#pkg"activate ."
using FourierGPE
gr(titlefontsize=12,size=(500,300),transpose=true,colorbar=false)

# ==== Units: ========================
# this example works in oscillator units

# ==== define user parameters =======
@with_kw mutable struct Params <: UserParams @deftype Float64
    # user parameters:
    κ = 0.1
end
par = Params()

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

# turn off the potential
V(x,y,t) = 0.0

#TODO define scaling
λx(t) = t
λy(t) = t

#TODO implement 2D scaling GPE
# scale espec (broadcast)
# scale interactions

import FourierGPE.Lgp!, FourierGPE.nlin!

function nlin!(dϕ,ϕ,sim::Sim{2},t)
    @unpack g,X = sim; x,y = X
    dϕ .= ϕ
    xspace!(dϕ,sim)
    @. dϕ *= g*abs2(dϕ) + V(x,y',t)
    kspace!(dϕ,sim)
    return nothing
end

function Lgp!(dϕ,ϕ,sim,t)
    @unpack γ,μ,espec = sim
    nlin!(dϕ,ϕ,sim,t)
    @. dϕ = -im*(1.0 - im*γ)*(dϕ + (espec - μ)*ϕ)
    return nothing
end


# ==== set simulation parameters ====
γ = 0.0
tf = 10.0
t = LinRange(ti,tf,Nt)
ϕi = ϕg

pare = Params()
sime = Sim(L,N,pare)
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

gif(anim,"./examples/scaleexpand.gif",fps=30)