using Test, SpecialFunctions, VortexDistributions
using LazyArrays, FillArrays
using Revise, FourierGPE


#--- Initialize simulation
# harmonic oscillator units
L = (18.0,18.0)
N = (256,256)
sim = Sim(L,N)
@unpack_Sim sim

# set simulation parameters
μ = 20.0

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
R(w) = sqrt(2*μ/w^2)
R(1)
rv = .8
healinglength(x,y,μ,g) = 1/sqrt(g*abs2(ψ0(x,y,μ,g)))
ξ0 = healinglength.(0.,0.,μ,g)
ξ = healinglength(rv,0.,μ,g)

pv = PointVortex(rv,0.,1)
nv = PointVortex(-rv,0,-1)

# dipole is a bit better behaved for testing
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

# Npad = N[1]/2 |> Int
# z0 = Zeros(ψi[:,1:Npad])
# ψc = Hcat(z0,view(ψi,:,:),z0)
# z0 = Zeros(ψc[1:Npad,:])
# ψc = Vcat(z0,ψc,z0) |> Matrix
# ϕc = fft(ψc)
# A = ifft(abs2.(ϕc)) |> fftshift
#
# heatmap(abs2.(ψi))
# heatmap(abs.(ϕi))
# heatmap(abs.(A))

#--- method to pad 2d array to twice size with zeros.

function zeropad(a)
    s = size(a)
    (isodd.(s) |> any) && error("Array dims must be divisible by 2")
    S = @. 2 * s
    t = @. s / 2 |> Int
    z = Zeros{eltype(a)}(t...)
    M = Hcat([z; z], a, [z; z])
    return Vcat([z z z z], M, [z z z z]) |> Matrix
end

function log10range(a,b,n)
    x = LinRange(log10(a),log10(b),n)
    return @. 10^x
end

function kespectrum(kp,ψ,x,y)
    dx,dy = diff(x)[1],diff(y)[1]
    Nx = 2*length(x)
    Lx = x[end]-x[1] + dx
    xp = LinRange(-Lx,Lx,Nx+1)[1:Nx]
    yp = xp
    ρ = @. sqrt(xp^2 + yp'^2)
    ψc = zeropad(ψ)
    ϕc = fft(ψc)
    A = ifft(abs2.(ϕc)) |> fftshift

    Ek = zero(kp)
    for i in eachindex(kp)
        k = kp[i]
        Ek[i]  = 0.5*k^3*sum(@. besselj0(k*ρ)*A)*dx*dy |> real
    end
    return Ek
end

ψc = zeropad(ψi)
ϕc = fft(ψc)
A = ifft(abs2.(ϕc)) |> fftshift

kmin = 0.1 #0.5*2*pi/R(1)
kmax = 2*pi/ξ0
Np = 200
kp = log10range(kmin,kmax,Np)

Ek = kespectrum(kp,ψi,x,y)
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
Wix = zeropad(wix)
Wiy = zeropad(wiy)
wixk = fft(Wix)
wiyk = fft(Wiy)
Cix = ifft(abs2.(wixk)) |> fftshift
Ciy = ifft(abs2.(wiyk)) |> fftshift
Ci = Cix .+ Ciy
heatmap(abs.(Ci))



d = 2*rv
kR = 2*pi/R(1)
kξ = 2*pi/ξ
kd = 2*pi/d
ka = 2*pi #k value associated with oscillator length
kmin = 0.1*kR
kmax = 2kξ
Np = 300
kp = log10range(kmin,kmax,Np)

function ikspectrum(kp,wconv,x,y)
    dx,dy = diff(x)[1],diff(y)[1]
    Nx = 2*length(x)
    Lx = x[end]-x[1] + dx
    xp = LinRange(-Lx,Lx,Nx+1)[1:Nx]
    yp = xp
    ρ = @. sqrt(xp^2 + yp'^2)

    Eki = zero(kp)
    for i in eachindex(kp)
        k = kp[i]
        Eki[i]  = 0.5*k*sum(@. besselj0(k*ρ)*wconv)*dx*dy |> real
    end
    return Eki
end

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
