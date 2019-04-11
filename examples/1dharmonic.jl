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
L = (20.0,)
N = (128,)
X,K,dX,dK,DX,DK,T = maketransforms(L,N)
espec = 0.5*k2(L...,N...)

# ====== Initialize simulation ======
sim = Sim(L,N,par)
@pack! sim = T,X,K,espec
initsim!(sim)
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
ψi .+= (randn(N...) |> complex)
ϕi = kspace(ψi,sim)
@pack! sim = ϕi
sim

# ====== Evolve in k space ==========
sol = runsim(sim.ϕi,sim)
# ===================================

# pull out the ground state:
ϕg = sol[end]
ψg = xspace(ϕg,sim)
showpsi(x,ψg)

# Add soliton


# ==== set simulation parameters ====
γ = 0.0
tf = Tv
t = LinRange(ti,tf,Nt)
ϕi = kspace(ψv,sim)
@pack! sim = tf,t,γ,ϕi
initsim!(sim)

@unpack_Sim sim
# ===================================

# ====== Evolve in k space ==========
solv = runsim(sim.ϕi,sim)
# ===================================

ϕf = solv[200]
ψf = xspace(ϕf,sim)
showpsi(ψf,x,y)

anim = @animate for i=1:Nt
    ψ = xspace(solv[i],sim)
    showpsi(x,ψ)
end

gif(anim,"./examples/soliton.gif",fps=30)
