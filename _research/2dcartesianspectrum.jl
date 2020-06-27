##
using Test, SpecialFunctions, VortexDistributions
using FourierGPE

## test in VSCode
using Pkg;Pkg.activate(".")
Pkg.test()

## simple diffeq test
using DifferentialEquations
f(u,p,t) = 1.01*u
u0 = 1/2
tspan = (0.0,1.0)
prob = ODEProblem(f,u0,tspan)
sol = solve(prob, Tsit5(), reltol=1e-8, abstol=1e-8)

## Initialize simulation
# harmonic oscillator units
L = (18.0,18.0)
N = (256,256)
sim = Sim(L,N)
@unpack_Sim sim

# set simulation parameters
μ = 20.0

## Time dependent potential function (here trivial t dep)
import FourierGPE.V
V(x,y,t) = 0.5*(x^2 + y^2)
ψ0(x,y,μ,g) = sqrt(μ/g)*sqrt(max(1.0-V(x,y,0.0)/μ,0.0)+im*0.0)

## make initial state
x,y = X
ψi = ψ0.(x,y',μ,g)
ϕi = kspace(ψi,sim)
@pack_Sim! sim
dx,dy = diff(x)[1],diff(y)[1]

## evolve
sol = runsim(sim)

## pull out ground state
ϕg = sol[end]
ψg = xspace(ϕg,sim)
showpsi(x,y,ψg)

## make dipole
R(w) = sqrt(2*μ/w^2)
R(1)
rv = .8
healinglength(x,y,μ,g) = 1/sqrt(g*abs2(ψ0(x,y,μ,g)))
ξ0 = healinglength.(0.,0.,μ,g)
ξ = healinglength(rv,0.,μ,g)

pv = PointVortex(rv,0.,1)
nv = PointVortex(-rv,0,-1)
v1 = ScalarVortex(ξ,pv)
v2 = ScalarVortex(ξ,nv)

psi = Torus(copy(ψg),x,y)
vortex!(psi,[v1;v2])
showpsi(x,y,psi.ψ)

## kespectrum
ψi .= psi.ψ
ϕi = kspace(ψi,sim) |> fftshift
kx,ky = K .|> fftshift
x,y = X
 
kmin = 0.1 #0.5*2*pi/R(1)
kmax = 2*pi/ξ0
Np = 200
kp = log10range(kmin,kmax,Np)

Ek = kespectrum(kp,ψi,x,y)
plot(kp,Ek,scale=:log10)



d = 2*rv
kR = 2*pi/R(1)
kξ = 2*pi/ξ
kd = 2*pi/d
ka = 2*pi #k value associated with oscillator length
kmin = 0.1*kR
kmax = 2kξ
Np = 300
kp = log10range(kmin,kmax,Np)

Eki = ikspectrum(kp,Ci,x,y)

kxi = kp*ka/kξ
plot(kxi,Eki,scale=:log10,grid=false,label=L"E_i(k)",legend=:topright)
plot!(kxi,2e7*kxi.^(-3),label=L"k^{-3}")
plot!(kxi,2e6*kxi,label=L"k")
ylims!(10^2,10^8)
xlims!(0.02,10)
vline!([kR*ka/kξ],ls=:dash,label=L"k_R")
vline!([kξ*ka/kξ],ls=:dash,label=L"k_\xi")
vline!([kd*ka/kξ],ls=:dash,label=L"k_d")
xlabel!(L"k\xi")
ylabel!(L"E_i(k)")
