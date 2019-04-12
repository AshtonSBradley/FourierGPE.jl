using Plots, LaTeXStrings, Pkg, Revise
gr(legend=false,titlefontsize=12,size=(500,300),transpose=true,colorbar=false)

#pkg"activate ."
using FourierGPE

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
    # user parameters:
    κ = 0.1
end
par = Params()

# ==== set simulation parameters ====
L = (60.0,)
N = (512,)
μ = 25.0
# old init code
# X,K,dX,dK,DX,DK,T = maketransforms(L,N)
# espec = 0.5*k2(L...,N...)
# @pack! sim = T,X,K,espec
# initsim!(sim)
# ====== Initialize simulation ======
sim = Sim(L,N,par)
@pack! sim = μ
@unpack_Sim sim
# ===================================

# declare the potential function
import FourierGPE.V
V(x,t) = 0.5*x^2
V(x,y,t)::Float64 = 0.5*(x^2 + y^2)
V(x,y,z,t)::Float64 = 0.5*(x^2 + y^2 + z^2)

# useful TF state
ψ0(x,μ,g) = sqrt(μ/g)*sqrt(max(1.0-V(x,0.0)/μ,0.0)+im*0.0)

x = X[1]
#make initial state
ψi = ψ0.(x,μ,g)
ϕi = kspace(ψi,sim)
@pack! sim = ϕi
sim

# ====== Evolve in k space ==========
sol = runsim(sim.ϕi,sim)
# ===================================

# pull out the ground state:
ϕg = sol[end]
ψg = xspace(ϕg,sim)
plot(x,g*abs2.(ψg),fill=(0,0.2),size=(600,200))
plot!(x,one.(x)*μ)
plot!(x,V.(x,0.0))
xlims!(-10,10); ylims!(0,1.3*μ)
title!(L"\textrm{local}\; \mu(x)")
xlabel!(L"x/a_x"); ylabel!(L"\mu(x)/\hbar\omega_x")

# Imprint dark soliton
ψf = xspace(sol[end],sim)
c = sqrt(μ)
ξ = 1/c
v = 0.5*c
xs = 0.
f = sqrt(1-(v/c)^2)

ψs = ψf.*(f*tanh.(f*(x .-xs)/ξ).+im*v/c);
showpsi(x,ψs)
xlims!(-10,10)
# ==== set simulation parameters ====
γ = 0.0
tf = 8*pi/sqrt(2); t = LinRange(ti,tf,Nt)
dt = 0.01π/μ
simSoliton = Sim(sim;γ=γ,tf=tf,t=t)
ϕi = kspace(ψs,simSoliton)
@pack! simSoliton = ϕi

@unpack_Sim simSoliton
# ===================================

# ====== Evolve in k space ==========
sols = runsim(simSoliton.ϕi,simSoliton)
# ===================================

ϕf = sols[end-4]
ψf = xspace(ϕf,simSoliton)
showpsi(x,ψf)

anim = @animate for i in 1:length(t)-4
    ψ = xspace(sols[i],simSoliton)
    y = g*abs2.(ψ)
    plot(x,y,fill=(0,0.2),size=(600,300),grid=false)
    xlims!(-10,10); ylims!(0,1.3*μ)
    title!(L"\textrm{local}\; \mu(x)")
    xlabel!(L"x/a_x"); ylabel!(L"\mu(x)/\hbar\omega_x")
end

gif(anim,"./examples/soliton.gif",fps=30)
