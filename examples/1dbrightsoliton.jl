using Plots, LaTeXStrings, Pkg, Revise
gr(legend=false,titlefontsize=12,size=(500,300),transpose=true,colorbar=false)

#pkg"activate ."
using FourierGPE

# declare the potential function
import FourierGPE.V
V(x,t) = zero(x) |> complex

# ==== Units: ========================
# this example works in oscillator units
function showpsi(x,ψ)
    p1 = plot(x,abs2.(ψ))
    xlabel!(L"x/a_x");ylabel!(L"|\psi|^2")
    p2 = plot(x,angle.(ψ))
    xlabel!(L"x/a_x");ylabel!(L"\textrm{phase}(\psi)")
    p = plot(p1,p2,layout=(2,1),size=(600,400))
    return p
end

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
us = 30
tf = π |> Float64
Nt = 400
t = LinRange(ti,tf,Nt)
# ====== Initialize simulation ======
sim = Sim(L,N,par)
@pack! sim = μ,g,γ,t,tf,Nt
@unpack_Sim sim
# ===================================
#soliton wavefunction
X,K = makearrays(L,N)
x = X[1]
ψs(x) = sqrt(Ns/2ξs)*sech(x/ξs)*exp(-im*us*x)
ψi = ψs.(x)
ϕi = kspace(ψi,sim)
@pack! sim = ϕi
sim
# ====== Evolve in k space ==========
sol = runsim(sim)
# ===================================

y = abs2.(sol[Nt-76])
plot(x,y,fill=(0, 0.2),size=(600,150),legend=false,grid=false,xticks=false,yticks=false,axis=false)

anim = @animate for i=1:2:Nt-76
    ψ = xspace(sol[i],sim)
    y=abs2.(ψ)
    plot(x,y,fill=(0, 0.2),size=(600,150),legend=false,grid=false,xticks=false,yticks=false,axis=false)
end

gif(anim, "./examples/brightsoliton.gif", fps = 24)

# anim = @animate for i in 1:length(t)-4
#     ψ = xspace(sols[i],simSoliton)
#     y = g*abs2.(ψ)
#     plot(x,y,fill=(0,0.2),size=(600,300))
#     xlims!(-10,10); ylims!(0,1.3*μ)
#     title!(L"\textrm{local}\; \mu(x)")
#     xlabel!(L"x/a_x"); ylabel!(L"\mu(x)/\hbar\omega_x")
# end
#
# gif(anim,"./examples/soliton.gif",fps=30)
