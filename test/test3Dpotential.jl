using Test, Parameters, BenchmarkTools

using FourierGPE

import FourierGPE:V
# ==== Units: ========================
# Units of ξ for length, 1/μ for time
# related by μ ≡ ħ²/mξ² = mc^2 where
# c = ħ/mξ is the speed of sound of
# the uniform system.

macro staticpotential(sim,V0)
    @eval sim = Sim(sim,V0=StaticPotential($V0))
end

# macro potential(sim,Vex)
#     @eval V(x,y,z,t)=$Vex
#     # @eval sim3 = Sim($sim,V=Potential($esc(ex)))
# end

# ==== set simulation parameters ====
L = (16.,16.,16.)
N = (64,64,64)
γ = 0.05
tf = 4/γ
Nt = 200
t = LinRange(0.,tf,Nt)
# ========= Initialize simulation ======
X = xvecs(L,N)
x,y,z = X
y = y'
z = reshape(z,1,1,length(z))

sim = Sim(L,N)
@staticpotential sim zero(x.*y.*z)

# #TODO get the macro to pass V to sim
# ex = :(x^11)
# dump(ex)
# @potential sim x^11
# sim = Sim(sim,V=Potential(ex))
# V(.1,0,0,0)
# Vtest(x,y,z,t)=sim.V.V
# Vtest(0.1,0,0,0)

V(x,y,z,t) = 0.25*(x^2 + y^2 + z^2)*(1+0.01*sin(t))

@pack! sim = γ,tf,Nt,t
@unpack μ,g = sim

# ========= useful state functions

ψ0(x,y,z,μ,g) = sqrt(μ/g)*sqrt(max(1.0-V(x,y,z,0)/μ,0.0)+im*0.0)
healing(x,y,z,μ,g) = 1/sqrt(g*abs2(ψ0(x,y,z,μ,g)))

#make initial state
ψi = ψ0.(x,y,reshape(z,1,1,length(z)),μ,g)

dψi = similar(ψi)

ϕi = kspace(ψi,sim)
dϕi = similar(ϕi)
@pack! sim = ϕi,γ,tf,t

@time nlin!(dψi,ψi,sim,0.1)

# ================
# define static potential array

testStatic = StaticPotential(zero(x.*y.*z))
