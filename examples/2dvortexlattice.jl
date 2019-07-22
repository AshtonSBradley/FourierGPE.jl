using Pkg, Revise, FourierGPE
gr(titlefontsize=12,size=(500,300),transpose=true,colorbar=false)

# ==== set simulation parameters ====
L = (25.0,25.0)
N = (256,256)
μ = 25.0

# ========= Initialize simulation ======
sim = Sim(L,N)
@pack! sim = μ
@unpack_Sim sim


# ===================================
# Two ways to set potential:

# Time dependent potential function (here trivial t dep)
import FourierGPE.V
V(x,y,t) = 0.5*(x^2 + y^2)
ψ0(x,y,μ,g) = sqrt(μ/g)*sqrt(max(1.0-V(x,y,0.0)/μ,0.0)+im*0.0)

# As a static Potential
# x,y = X
# Vs(x,y) = 0.5*(x^2 + y^2)
# ψ0(x,y,μ,g) = sqrt(μ/g)*sqrt(max(1.0-Vs(x,y)/μ,0.0)+im*0.0)
# sim = Sim(sim,V0=Vs.(x,y'))

# ==== make initial state
x,y = X
ψi = ψ0.(x,y',μ,g)
ϕi = kspace(ψi,sim)
@pack! sim = ϕi,μ

# ==== Evolve in k space
@time sol = runsim(sim)

# ==== ground state
ϕg = sol[end]
ψg = xspace(ϕg,sim)
showpsi(x,y,ψg)

# ==== set simulation parameters
γ = 0.3
t = LinRange(ti,tf,Nt)
ϕi = kspace(ψg,sim)
# reltol = 1e-7
# alg = DP5()
@pack! sim = tf,t,γ,ϕi #,reltol,alg
@unpack_Sim sim

import FourierGPE: nlin!, Lgp!
sim

function Lz!(dϕ,ϕ,sim::Sim{2},t)
    #TODO
    @unpack X,K,T,params = sim; x,y = X; kx, ky = K
    χ = T.psi.x[1]
    Ω = params.κ

    # x*ky
    χ .= ϕ
    @. χ *= ky'
    xspace!(χ,sim)
    @. χ *= x
    kspace!(χ,sim)
    @. dϕ += Ω*χ

    # -y*kx
    χ .= ϕ
    @. χ *= kx
    xspace!(χ,sim)
    @. χ *= y'
    kspace!(χ,sim)
    @. dϕ -= Ω*χ
    return nothing
end

function Lz(ϕ,sim::Sim{2})
    @unpack X,K,params = sim; x,y = X; kx, ky = K
    Ω = params.κ
    # x*ky
    psi1 = ϕ.*ky'
    xspace!(psi1,sim)
    @. psi1 *= Ω*x
    kspace!(psi1,sim)

    # y*kx
    psi2 = ϕ.*kx
    xspace!(psi2,sim)
    @. psi2 *= Ω*y'
    kspace!(psi2,sim)
    return psi1 .- psi2
end

function nlin!(dϕ,ϕ,sim::Sim{2},t)
    @unpack g,X,V0 = sim; x,y = X
    y = y'
    dϕ .= ϕ
    xspace!(dϕ,sim)
    @. dϕ *= V0 + V(x,y,t) + g*abs2(dϕ)
    kspace!(dϕ,sim)
    return nothing
end

function Lgp!(dϕ,ϕ,sim,t)
    @unpack γ,μ,espec = sim
    nlin!(dϕ,ϕ,sim,t)
    dϕ .+= Lz(ϕ,sim)
    @. dϕ = -im*(1.0 - im*γ)*(dϕ + (espec - μ)*ϕ)
    return nothing
end

params = Params(κ = 0.9)
@pack! sim = params
sim
ti = 0.; tf = 5.
t = LinRange(ti,tf,Nt)
ϕi = ϕi + randn(N)+im*randn(N)
# ==== Evolve in k space
solv = runsim(sim)

# ===================================
ϕf = solv[end]
ψf = xspace(ϕf,sim)
showpsi(x,y,ψf)


anim = @animate for i in eachindex(t)
    ψ = xspace(solv[i],sim)
    showpsi(x,y,ψ)
end

gif(anim,"./examples/vortexlattice.gif",fps=25)
