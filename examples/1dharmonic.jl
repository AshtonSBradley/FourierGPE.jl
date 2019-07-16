using Plots, LaTeXStrings, Pkg, Revise
gr(grid=false,legend=false,titlefontsize=12,size=(500,300),transpose=true,colorbar=false)

using FourierGPE

# ==== set simulation parameters ====
L = (40.0,)
N = (512,)
μ = 25.0

sim = Sim(L,N)
@pack! sim = μ
@unpack_Sim sim

# ===================================

# declare the potential function
import FourierGPE.V
V(x,t) = 0.5*x^2

# useful TF state
ψ0(x,μ,g) = sqrt(μ/g)*sqrt(max(1.0-V(x,0.0)/μ,0.0)+im*0.0)

#make initial state
x = X[1]
ψi = ψ0.(x,μ,g)
ϕi = kspace(ψi,sim)
@pack_Sim! sim

# ====== Evolve in k space ==========
sol = runsim(sim)
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
simSoliton = Sim(sim;γ=γ,tf=tf,t=t)
ϕi = kspace(ψs,simSoliton)
@pack! simSoliton = ϕi

#update all variables in current global workspace
@unpack_Sim simSoliton
# ===================================

# ====== Evolve in k space ==========
sols = runsim(simSoliton)
# ===================================

ϕf = sols[end-3]
ψf = xspace(ϕf,simSoliton)
showpsi(x,ψf)

anim = @animate for i in 1:length(t)-3
    ψ = xspace(sols[i],simSoliton)
    y = g*abs2.(ψ)
    plot(x,y,fill=(0,0.2),size=(600,300),grid=false)
    xlims!(-10,10); ylims!(0,1.3*μ)
    title!(L"\textrm{local}\; \mu(x)")
    xlabel!(L"x/a_x"); ylabel!(L"\mu(x)/\hbar\omega_x")
end

gif(anim,"./examples/soliton.gif",fps=30)
