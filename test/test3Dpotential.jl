using Test, Parameters, BenchmarkTools, Revise

using FourierGPE

# ==== Units: ========================
# Units of ξ for length, 1/μ for time
# related by μ ≡ ħ²/mξ² = mc^2 where
# c = ħ/mξ is the speed of sound of
# the uniform system.

macro potential(sim,V)
    fn = Symbol("V(x,y,z,t)")
    quote
        function $(esc(fn))()
            return $V
        end
        @eval sim = Sim(sim,V=Potential($V))
    end
end

macro staticpotential(sim,V0)
    @eval sim = Sim(sim,V0=StaticPotential($V0))
end

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

# par = Params()
# p1 = Potential(:(0.25*(x^2 + y^2 + z^2)*(1+t)))
# V0 = StaticPotential(zero(x.*y.*z))

sim = Sim(L,N)
@staticpotential sim zero(x.*y.*z)

#TODO get the macro to define the function
@potential sim :(0.25*(x^2 + y^2 + z^2)*(1+0.01*sin(t)))
V(x,y,z,t) = 0.25*(x^2 + y^2 + z^2)*(1+0.01*sin(t))

@pack! sim = γ,tf,Nt,t

@unpack μ,g = sim
# @unpack_Sim sim

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



@potential sim :(0.25*(x^2 + y^2 + 4*z^2))
@staticpotential sim zero(x.*y.*z)
