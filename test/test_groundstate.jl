# ground state tests

import FourierGPE:V
V(x,t) = 0.5*x^2


# ==== set simulation parameters ====
L = (40.0,)
N = (512,)
μ = 25.0

# === collect parameters ===========
sim = Sim(L,N)
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

psitest = xspace(sol.u[end],sim)
npeak = abs2.(psitest[256])
@test isapprox(g*npeak,μ,rtol=1e-2)
