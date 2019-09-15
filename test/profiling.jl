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
    # @. dϕ *= g*abs2(dϕ)
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
# factor of 10 speed up in MKL
# if there isn't a speed up, there must be another bottleneck

@btime kspace!(dϕi,sim)
@btime xspace!(dϕi,sim)
@btime @. dϕi *= g*abs2(dϕ)



y = abs2.(xspace(sol[Nt-76],sim))
plot(x,y,fill=(0, 0.2))

anim = @animate for i=1:Nt-8
    ψ = xspace(sol[i],sim)
    y = abs2.(ψ)
    plot(x,y,c=msg,fill=(0, 0.4,msg))
end

gif(anim, "./examples/brightsoliton.gif", fps = 25)
