using Plots, LaTeXStrings, Pkg, Revise
gr(colorbar=false,size=(600,150),legend=false,grid=false,xticks=true,yticks=true,axis=true)

using FourierGPE

# ==== set simulation parameters ====
L = (60.0,)
N = (512,)
μ = 25.0
g = 0.01
γ = 0.0
tf = 2.0pi
Nt = 150
ti = 0.0
t = LinRange(ti,tf,Nt)
# ====== Initialize simulation ======
sim = Sim(L,N)
@pack! sim = μ,g,γ,t,tf,Nt
@unpack_Sim sim
# ===================================
#bogoliubov state
X,K = makearrays(L,N)
x = X[1]
k = K[1]
# define Bog state

f(k) = 1 + 4*μ/k^2
u(k) = 0.5*(sqrt(f(k))+1)/f(k)^(1/4)
v(k) = 0.5*(sqrt(f(k))-1)/f(k)^(1/4)
lam = 0.000001
bog(x,k) = u(k)*exp(im*k*x) - conj(v(k))*exp(-im*k*x)
ψb(x,k) = sqrt(μ/g)*(1 + lam*bog(x,k))

kb = k[20]
ψi = ψb.(x,kb)
ϕi = kspace(ψi,sim)

# Set time evolution and pack
reltol = 1e-7
alg = Vern7()

@pack_Sim! sim

# ====== Evolve in k space ==========
sol = runsim(sim)
# ===================================

y = abs2.(xspace(sol[end],sim));plot(x,y)

#make a movie?
anim = @animate for i=1:Nt
    ψ = xspace(sol[i],sim)
    y=abs2.(ψ)
    plot(x,y)
end

gif(anim, "./examples/bogoliubov.gif", fps = 25)
