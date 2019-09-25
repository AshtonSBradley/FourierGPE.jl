using Test, SpecialFunctions, VortexDistributions
using LazyArrays, FillArrays
using Revise, FourierGPE

using Plots, ColorSchemes
c1 = cgrad(ColorSchemes.linear_blue_5_95_c73_n256.colors)
c2 = cgrad(ColorSchemes.turbo.colors)

import FourierGPE: showpsi
function showpsi(x, y, ψ)
    p1 = heatmap(
        x,
        y,
        abs2.(ψ),
        aspectratio = 1,
        c = c1,
        titlefontsize = 12,
        transpose = true,
        colorbar = false,
    )
    xlims!(x[1], x[end])
    ylims!(y[1], y[end])
    xlabel!(L"x")
    ylabel!(L"y")
    title!(L"|\psi|^2")
    p2 = heatmap(
        x,
        y,
        angle.(ψ),
        aspectratio = 1,
        c = c2,
        titlefontsize = 12,
        transpose = true,
        colorbar = false,
    )
    xlims!(x[1], x[end])
    ylims!(y[1], y[end])
    xlabel!(L"x")
    ylabel!(L"y")
    title!(L"\textrm{phase} (\psi)")
    p = plot(p1, p2, size = (600, 300))
    return p
end

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

kmin = 0.1 #0.5*2*pi/R(1)
kmax = 2*pi/ξ0
Np = 200
ks = LinRange(kmin,kmax,Np)
kp = @. log(exp(ks))

# Ek = zero(kp)
# Nx = 2*N[1]
# xp = LinRange(-L[1],L[1],Nx+1)[1:Nx]
# yp = xp
# ρ = @. sqrt(xp^2 + yp'^2)

function kspectrum(kp,ψ,x,y)
    dx,dy = diff(x)[1],diff(y)[1]
    N = length(x)
    L = last(x)-first(x) + dx
    Npad = N/2 |> Int
    z0 = Zeros(ψ[:,1:Npad])
    ψc = Hcat(z0,view(ψ,:,:),z0)
    z0 = Zeros(ψc[1:Npad,:])
    ψc = Vcat(z0,ψc,z0) |> Matrix
    ϕc = fft(ψc)
    A = ifft(abs2.(ϕc)) |> fftshift

    Ek = zero(kp)
    for i in eachindex(kp)
        k = kp[i]
        Nx = 2*N
        xp = LinRange(-L,L,Nx+1)[1:Nx]
        yp = xp
        ρ = @. sqrt(xp^2 + yp'^2)
        Ek[i] = 0.5*k^3*sum(@. besselj0(k*ρ)*A)*dx*dy |> real
    end
    return Ek
end

Ek = kspectrum(kp,ψi,x,y)
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

function log10range(a,b,n)
    x = LinRange(log10(a),log10(b),n)
    return @. 10^x
end

d = 2*rv
kR = 2*pi/R(1)
kξ = 2*pi/ξ
kd = 2*pi/d
ka = 2*pi #k value associated with oscillator length
kmin = 0.1*kR
kmax = 2kξ
Np = 300
kp = log10range(kmin,kmax,Np)

Ei = zero(kp)
Nx = 2*N[1]
xp = LinRange(-L[1],L[1],Nx+1)[1:Nx]
yp = xp
ρ = @. sqrt(xp^2 + yp'^2)

for i in eachindex(kp)
    k = kp[i]
    Ei[i] = 0.5*k*sum(@. besselj0(k*ρ)*Ci)*dx*dy |> real
end

kxi = kp*ka/kξ
plot(kxi,Ei,scale=:log10,grid=false,label=L"E_i(k)",legend=:topright)
plot!(kxi,2e7*kxi.^(-3),label=L"k^{-3}")
plot!(kxi,2e6*kxi,label=L"k")
ylims!(10^2,10^8)
xlims!(0.02,10)
vline!([kR*ka/kξ],ls=:dash,label=L"k_R")
vline!([kξ*ka/kξ],ls=:dash,label=L"k_\xi")
vline!([kd*ka/kξ],ls=:dash,label=L"k_d")
xlabel!(L"k\xi")
ylabel!(L"E_i(k)")
