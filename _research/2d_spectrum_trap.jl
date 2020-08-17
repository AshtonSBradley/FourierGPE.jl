using VortexDistributions, FourierGPE, Plots, LaTeXStrings
using SpecialFunctions


## Initialize simulation
L = (18.0,18.0)
N = (256,256)
sim = Sim(L,N)
@unpack_Sim sim

## set simulation parameters
μ = 12.0
tf = 10.0
t = LinRange(ti,tf,Nt)

# Time dependent potential function (here trivial t dep)
import FourierGPE.V
V(x,y,t) = 0.5*(x^2 + y^2)
ψ0(x,y,μ,g) = sqrt(μ/g)*sqrt(max(1.0-V(x,y,0.0)/μ,0.0)+im*0.0)

## make initial state
x,y = X
ψi = ψ0.(x,y',μ,g)
ϕi = kspace(ψi,sim)
@pack_Sim! sim

## evolve
@time sol = runsim(sim)

## plot ground state
ϕg = sol[end]
ψg = xspace(ϕg,sim)
healinglength(x,y,μ,g) = 1/sqrt(g*abs2(ψ0(x,y,μ,g)))
R(μ,g) = sqrt(2*μ)
showpsi(x,y,ψg)
Rtf = R(μ,g)
n0 = maximum(abs2.(ψg))


## spectra of ground state 
ξ = healinglength(0,0,μ,g)
kR = 2*pi/Rtf
kξ = 2*pi

kmin = 0.1*kR
kmax = kξ

Np = 300
k = log10range(kmin,kmax,Np)

@time Ek = kespectrum(k,ψg,X,K)
plot(k*ξ,Ek,scale=:log10)

# analytic
f(x) = (sin(x)-x*cos(x))^2/x^3
Ektf(k,R) = π*n0*R*f.(k*R)

# linear scale
plot(k*ξ,Ek)
plot!(k*ξ,Ektf.(k,Rtf), c=:red)

## spectrum
kR = 2*pi/Rtf
kξ = 2*pi

kmin = 0.1*kR
kmax = 1.3kξ
Np = 200
k = log10range(kmin,kmax,Np)

Eki = ikespectrum(k,ψg,X,K)
Ekc = ckespectrum(k,ψg,X,K)
Ekq = qpespectrum(k,ψg,X,K)
plot(k*ξ,Eki,label=L"E_k^i(k)")
plot!(k*ξ,Ekc,label=L"E_k^c(k)")
plot!(k*ξ,Ekq,label=L"E_k^q(k)")

plot(k*ξ,Ekq,scale=:log10)
vline!([kR*ξ],ls=:dash)
vline!([1.0],ls=:dash)

## also plot interaction and trap energies without vortex
using QuadGK
# trap energy
int1(k)=quadgk(x-> x^2*sqrt(1-x^2)*besselj0(k*x),0,1)[1]

# interaction energy
int2(k)=quadgk(x-> x*(1-x^2)*besselj0(k*x),0,1)[1]

plot(k,)
## make central vortex
ψd = copy(ψg)
ξv = healinglength(0,0.,μ,g)
# can use default ξ = 1
xp,yp = 0.,d/2
xn,yn = 0.,-d/2
pv = PointVortex(xp,yp,1)
nv = PointVortex(xn,yn,-1)
dipole = ScalarVortex(ξv,[nv;pv])

psi = Torus(ψd,x,y) # set field topology for VortexDistributions
vortex!(psi,dipole) # make dipole
showpsi(x,y,ψd)

## evolve
tf = 1.0
t = LinRange(ti,tf,Nt)
γ = 0.
ψi .= psi.ψ
ϕi = kspace(ψi,sim)

@pack_Sim! sim
@time sol = runsim(sim)


## plot final state
ϕf = sol[end]
ψf = xspace(ϕf,sim)
showpsi(x,y,ψf)

## analyse
kL = 2*pi/L[1]
kξ = 2*pi

kmin = 0.1*kL
kmax = kξ

Np = 200
k = log10range(kmin,kmax,Np)

@time Ek = kespectrum(k,ψf,X,K)
plot(k,Ek,scale=:log10)

# heatmap(log.(abs2.(A |> fftshift) .+ eps.()))

## incompressible spectrum
kL = 2*pi/L[1]
kξ = 2*pi
kd = 2*pi/d

kmin = 0.1*kL
kmax = 1.3kξ
Np = 400
k = log10range(kmin,kmax,Np)

@time Eki = ikespectrum(k,ψf,X,K)

## power-law plot in logspace
n0 = abs2.(ψf[1,1])
E0 = π*n0
p1 = plot(k,Eki/E0,scale=:log10,grid=false,label=L"E_i(k)",legend=:bottomleft)
# plot!(k,200k.^(-3),label=L"k^{-3}")
# plot!(k,2000*k,label=L"k")
# ylims!(2e-3,1e2)
# xlims!(0.02,10)
vline!([kL],ls=:dash,label=L"k_L")
vline!([kd],ls=:dash,label=L"k_d")
vline!([kξ],ls=:dash,label=L"k_\xi")
xlabel!(L"k\xi")
ylabel!(L"E_i(k)")

# test against analytic form
# use the numerical dipole distance
# vort = findvortices(Torus(ψf,x,y)) |> rawData
# d = vort[:,2] |> diff
d = 2
Λ = 0.8249
f(x) = x*(besselk(1,x)*besseli(0,x)-besselk(0,x)*besseli(1,x))
F(x,Λ) = f(x/(2Λ))^2/x
F(x) = F(x,Λ)
Ed(k,d) = 2*F(k)*(1-besselj0(k*d))
plot!(k,Ed.(k,[d]),label=L"{\cal E}_i^a(k)")
