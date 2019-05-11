# ground state tests
using Test, Revise, FourierGPE
@with_kw mutable struct Params <: UserParams @deftype Float64
    # user parameters: add with defaults.
    ω = 2*pi
    V::Expr = :(V(x,t) = 0.5*ω^2*x^2 )
end

# ==== set simulation parameters ====
L = (40.0,)
N = (512,)
μ = 25.0
ω = 1.0
Vdef = :(V(x,t) = 0.5*ω^2*x^2)

# === collect parameters ===========
eval(Vdef)
par = Params(ω = ω,V=Vdef)

sim = Sim(L,N,par)
@pack! sim = μ
@unpack_Sim sim
# ===================================

# useful TF state
ψ0(x,μ,g) = sqrt(μ/g)*sqrt(max(1.0-V(x,0.0)/μ,0.0)+im*0.0)

#evaluate on grid
x = X[1]
ψi = ψ0.(x,μ,g)

#make initial state
ϕi = kspace(ψi,sim)
@pack! sim = ϕi

# ====== Evolve in k space without error ==========

runsim(sim;info=false)

sol,err = testsim(sim)
@test err == false

psitest = xspace(sol[end],sim)
npeak = abs2.(psitest[256])
@test isapprox(g*npeak,μ,rtol=1e-2)
