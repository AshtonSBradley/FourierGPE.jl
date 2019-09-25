using DSP, Test, SpecialFunctions, VortexDistributions
using Revise, FourierGPE

#--- Initialize simulation
L = (18.0,18.0)
N = (256,256)
sim = Sim(L,N)
@unpack_Sim sim

# set simulation parameters
μ = 12.0

# Time dependent potential function (here trivial t dep)
import FourierGPE.V
V(x,y,t) = 0.5*(x^2 + y^2)
ψ0(x,y,μ,g) = sqrt(μ/g)*sqrt(max(1.0-V(x,y,0.0)/μ,0.0)+im*0.0)

# make initial state
x,y = X
ψi = ψ0.(x,y',μ,g)
ϕi = kspace(ψi,sim)
@pack_Sim! sim
dx,dy = diff(x)[1],diff(y)[1]
#---

#--- evolve
@time sol = runsim(sim)

# ground state

ϕg = sol[end]
ψg = xspace(ϕg,sim)
showpsi(x,y,ψg)
#---


#---
using VortexDistributions
R(w) = sqrt(2*μ/w^2)
R(1)
rv = 0.
healinglength(x,y,μ,g) = 1/sqrt(g*abs2(ψ0(x,y,μ,g)))
ξ0 = healinglength.(0.,0.,μ,g)
ξ = healinglength(rv,0.,μ,g)

pv = PointVortex(rv,0.,1)
vi = ScalarVortex(ξ,pv)


psi = Torus(copy(ψg),x,y)
vortex!(psi,vi)
showpsi(x,y,psi.ψ)

ψi .= psi.ψ
ϕi = kspace(ψi,sim) |> fftshift
kx,ky = K .|> fftshift
x,y = X
#---

#--- construct polar spectrum
# first convolve, then bessel
A = -conv(ψi,conj.(ψi))
heatmap(abs2.(ψi))
heatmap(abs.(ϕi))
heatmap(abs.(A))

kmin = 0.1;#0.5*2*pi/R(1)
kmax = 2*pi/ξ0/2
Np = 100
ks = LinRange(kmin,kmax,Np)
kp = @. log(exp(ks))

Ek = zero(kp)
Nx = 2*N[1]
xp = LinRange(-L[1],L[1],Nx)[1:Nx-1]
yp = xp
ρ = @. sqrt(xp^2 + yp'^2)

for i in eachindex(kp)
    k = kp[i]
    Ek[i] = 0.5*k^3*sum(@. besselj0(k*ρ)*A)*dx*dy |> real
end


plot(kp,Ek .+ eps.() .+ abs(minimum(Ek)),yscale=:log10,xscale=:log10)
#---

#--- test incompressible spectrum
k2field = k2(K)
psi = XField(ψi,X,K,k2field)
vx,vy = velocity(psi)
rho = abs2.(ψi)
wx = @. sqrt(rho)*vx; wy = @. sqrt(rho)*vy
Wi, Wc = helmholtz(wx,wy,psi)

# convolutions
Cix = -conv(conj.(Wi[1]),Wi[1])
Ciy = -conv(conj.(Wi[2]),Wi[2])
Ci = Cix .+ Ciy
heatmap(abs.(Ci))

Ei = zero(kp)
Nx = 2*N[1]
xp = LinRange(-L[1],L[1],Nx)[1:Nx-1]
yp = xp
ρ = @. sqrt(xp^2 + yp'^2)

for i in eachindex(kp)
    k = kp[i]
    Ei[i] = 0.5*k*sum(@. besselj0(k*ρ)*Ci)*dx*dy |> real
end
plot(kp,Ei)
