# ground state tests
using FourierGPE, Test
@with_kw mutable struct Params <: UserParams @deftype Float64
    # parameters (at least a placeholder):
    κ = 0.1
end
par = Params()

# ==== set simulation parameters ====
L = (40.0,)
N = (512,)
μ = 25.0

#X,K,dX,dK = makearrays(L,N)
#
# T = makeT(X,K)


sim = Sim(L,N,par)
@pack! sim = μ
@unpack_Sim sim
# ===================================

# declare the potential function
import FourierGPE.V
V(x,t) = 0.5*x^2

# useful TF state
ψ0(x,μ,g) = sqrt(μ/g)*sqrt(max(1.0-V(x,0.0)/μ,0.0)+im*0.0)

#evaluate on grid
x = X[1]
ψi = ψ0.(x,μ,g)

# using Plots
# plot(x,g*abs2.(ψi))

#make initial state
ϕi = kspace(ψi,sim)
@pack! sim = ϕi

# ====== Evolve in k space without error ==========

function testsim(sim)
    err = false
    sol = try
            runsim(sim;info=false)
        catch e
            err = true
        end
return sol,err
end

runsim(sim;info=false,nfiles=false)

sol,err = testsim(sim)
@test err == false

psitest = xspace(sol[end],sim)
npeak = abs2.(psitest[256])
@test isapprox(g*npeak,μ,rtol=1e-2)
