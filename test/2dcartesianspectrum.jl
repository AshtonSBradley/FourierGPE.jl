using DSP, Test, SpecialFunctions, VortexDistributions
using LazyArrays, FillArrays
using Revise, FourierGPE

#--- Initialize simulation
# harmonic oscillator units
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
rv = 1.5
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

ψi .= psi.ψ
ϕi = kspace(ψi,sim) |> fftshift
kx,ky = K .|> fftshift
x,y = X
#---

#--- construct polar spectrum
# convolve, then bessel

Npad = N[1]/2 |> Int
z0 = Zeros(ψi[:,1:Npad])
ψc = Hcat(z0,ψi,z0)
z0 = Zeros(ψc[1:Npad,:])
ψc = Vcat(z0,ψc,z0) |> Matrix
ϕc = fft(ψc)
A = ifft(abs2.(ϕc)) |> fftshift

heatmap(abs2.(ψi))
heatmap(abs.(ϕi))
heatmap(abs.(A))

kmin = 0.1 #0.5*2*pi/R(1)
kmax = 2*pi/ξ0
Np = 200
ks = LinRange(kmin,kmax,Np)
kp = @. log(exp(ks))

Ek = zero(kp)
Nx = 2*N[1]
xp = LinRange(-L[1],L[1],Nx+1)[1:Nx]
yp = xp
ρ = @. sqrt(xp^2 + yp'^2)

for i in eachindex(kp)
    k = kp[i]
    Ek[i] = 0.5*k^3*sum(@. besselj0(k*ρ)*A)*dx*dy |> real
end

plot(kp,Ek,scale=:log10)

#---





#--- test incompressible spectrum
k2field = k2(K)
psi = XField(ψi,X,K,k2field)
vx,vy = velocity(psi)
rho = abs2.(ψi)
wx = @. sqrt(rho)*vx; wy = @. sqrt(rho)*vy
Wi, Wc = helmholtz(wx,wy,psi)
wix,wiy = Wi

# convolutions
z0 = Zeros(wix[:,1:Npad])
Wix = Hcat(z0,wix,z0)
Wiy = Hcat(z0,wiy,z0)
z0 = Zeros(Wix[1:Npad,:])
Wix = Vcat(z0,Wix,z0) |> Matrix
Wiy = Vcat(z0,Wiy,z0) |> Matrix
Wixk = fft(Wix)
Wiyk = fft(Wiy)
Cix = ifft(abs2.(Wixk)) |> fftshift
Ciy = ifft(abs2.(Wiyk)) |> fftshift
Ci = Cix .+ Ciy
heatmap(abs.(Ci))


kR = 0.5*2*pi/R(1)
kξ = 2*pi
kmin = 0.02
kmax = 1.5*kξ
Np = 400
ks = LinRange(kmin,kmax,Np)
kp = @. log(exp(ks))
Ei = zero(kp)
Nx = 2*N[1]
xp = LinRange(-L[1],L[1],Nx+1)[1:Nx]
yp = xp
ρ = @. sqrt(xp^2 + yp'^2)

for i in eachindex(kp)
    k = kp[i]
    Ei[i] = 0.5*k*sum(@. besselj0(k*ρ)*Ci)*dx*dy |> real
end

plot(kp,Ei,scale=:log10,grid=false,label=L"E_i(k)",legend=:topright)
plot!(kp,2e7kp.^(-3),label=L"k^{-3}")
plot!(kp,2e6kp,label=L"k")
ylims!(10^2,10^7.5)
vline!([kR],ls=:dash,label=L"k_R")
vline!([kξ],ls=:dash,label=L"k_\xi")
xlabel!(L"k")
ylabel!(L"E_i(k)")
