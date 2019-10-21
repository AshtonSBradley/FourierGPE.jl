using Revise, FourierGPE
gr(titlefontsize=12,size=(500,300),transpose=true,colorbar=false)

# Units
# this example works in oscillator units

function jrphase(x,y)
    return 4x/(1.5 + (2y^2)+x^2)
end
#--- Initialize simulation

L = (30.0,30.0)
N = (128,128)
sim = Sim(L,N)
@unpack_Sim sim

#--- simulation parameters
μ = 60.0
ξ = 1/sqrt(μ)
offset = 2.0
# potential
import FourierGPE.V
V1 = 5*μ
sig = 3*ξ
V(x,y,t) = 0.5*(x^2 + y^2) + V1*exp(-((x-offset)^2 + y^2)/sig^2)

# TF state
ψ0(x,y,μ,g) = sqrt(μ/g)*sqrt(max(1.0-V(x,y,0.0)/μ,0.0)+im*0.0)

# make initial state
x,y = X
ψi = ψ0.(x,y',μ,g)
ϕi = kspace(ψi,sim)

@pack_Sim! sim

#--- evolve and plot final state
sol = runsim(sim)
ϕg = sol[end]
ψg = xspace(ϕg,sim)
showpsi(x,y,ψg)



#--- imprint and evolve JR for short time
# shift = sig
# ψjr = @. ψg*exp(-im*jrphase((x - offset-shift)/ξ,y'/ξ))
# showpsi(x,y,ψjr)
# ϕi = kspace(ψjr,sim)
# γ = 0.0
# tf = .3
# t = LinRange(ti,tf,Nt)
#
# simjr = Sim(L,N)
# @pack! simjr = tf,t,γ,ϕi

#--- try jump of obstacle
shift = -2*sig
# ψjr = @. ψg*exp(-im*jrphase((x - offset-shift)/ξ,y'/ξ))

V(x,y,t) = 0.5*(x^2 + y^2) + V1*exp(-((x-offset-shift)^2 + y^2)/sig^2)

ϕi = kspace(ψg,sim)
γ = 0.0
tf = .3
t = LinRange(ti,tf,Nt)

simjr = Sim(L,N)
@pack! simjr = tf,t,γ,ϕi

#--- evolve
soljr = runsim(simjr)

ϕf = soljr[end]
ψf = xspace(ϕf,simjr)
showpsi(x,y,ψf)

anim = @animate for i=1:Nt-6
    ψ = xspace(soljr[i],simjr)
    showpsi(x,y,ψ)
end

gif(anim,"./examples/jr.gif",fps=30)
