using Plots, LaTeXStrings, Pkg, Revise
gr(colorbar=false,size=(600,150),legend=false,grid=false,xticks=false,yticks=false,axis=false)

using FourierGPE

# declare the potential function
import FourierGPE.V
V(x,t) = zero(x) |> complex

# ==== define user parameters =======
@with_kw mutable struct Params <: UserParams @deftype Float64
    # parameters (at least a placeholder):
    κ = 0.1
end
par = Params()

# ==== set simulation parameters ====
L = (60.0,)
N = (512,)
μ = 25.0
g = -0.01
γ = 0.0
Ns = 200
ξs = 2/abs(g)/Ns
us = 20
tf = π |> Float64
Nt = 150
t = LinRange(0.,tf,Nt)

# ====== Initialize simulation ======
sim = Sim(L,N,par)
@pack! sim = μ,g,γ,t,tf,Nt
@unpack_Sim sim
# ===================================
#soliton wavefunction
X,K = makearrays(L,N)
x = X[1]
ψs(x) = sqrt(Ns/2ξs)*sech(x/ξs)*exp(im*us*x)
ψi = ψs.(x)
ϕi = kspace(ψi,sim)
@pack_Sim! sim
sim
# ====== Evolve in k space ==========
sol = runsim(sim)
# ===================================

y = abs2.(xspace(sol[Nt-76],sim))
plot(x,y,fill=(0, 0.2))

anim = @animate for i=1:Nt-8
    ψ = xspace(sol[i],sim)
    y=abs2.(ψ)
    plot(x,y,fill=(0, 0.2))
end

gif(anim, "./examples/brightsoliton.gif", fps = 25)
