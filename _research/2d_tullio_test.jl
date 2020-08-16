using VortexDistributions, FourierGPE, Plots, LaTeXStrings, SpecialFunctions
using Tullio, LoopVectorization

## Initialize simulation
# units of healing length, chemical potential
L = (100.,100.)
N = (512,512)
sim = Sim(L,N)
@unpack_Sim sim

## homogeneous state
μ = 1.0
g = 0.01
ψ0(x,y,μ,g) = sqrt(μ/g)*sqrt(max(1.0-V(x,y,0.0)/μ,0.0)+im*0.0)
healinglength(x,y,μ,g) = 1/sqrt(g*abs2(ψ0(x,y,μ,g)))
x,y = X
kx,ky = K
ψg = ψ0.(x,y',μ,g)
showpsi(x,y,ψg)

## initial dipole, with periodic phase
ψd = copy(ψg)
d = 12
ξv = healinglength(0.,0.,μ,g)
# can use default ξ = 1
xp,yp = 0.,d/2
xn,yn = 0.,-d/2
pv = PointVortex(xp,yp,1)
nv = PointVortex(xn,yn,-1)
dipole = [nv;pv]

psi = Torus(ψd,x,y) # set field topology for VortexDistributions
vortex!(psi,dipole) # make dipole

showpsi(x,y,ψd)

## test new methods
kd = 2*pi/d
kL = 2*pi/L[1]
kξ = 2*pi
kmin = 0.1*kL
kmax = kξ

Np = 300
k = log10range(kmin,kmax,Np)

import FourierGPE:bessel_reduce

function bessel_reduce(k,x,y,C)
    dx,dy = diff(x)[1],diff(y)[1]
    Nx = 2*length(x)
    Lx = x[end] - x[begin] + dx
    xp = LinRange(-Lx,Lx,Nx+1)[1:Nx]
    yp = xp
    ρ = hypot.(xp,yp')
    E = zero(k)
    @tullio E[i] = real(besselj0(k[i]*ρ[p,q])*C[p,q])
    @. E *= 0.5*k*dx*dy 
    return E 
end

## power-law plot in logspace
@time Eki = ikespectrum(k,ψd,X,K)

Eqp = qpespectrum(k,ψd,X,K)
Ekc = ckespectrum(k,ψd,X,K)

## plot
n0 = abs2.(ψd[1,1])
E0 = π*n0
p1 = plot(k,Eki/E0,scale=:log10,grid=false,label=L"\epsilon_I(k)",legend=:bottomleft)
plot!(k,Eqp/E0,label=L"\epsilon_Q(k)")
plot!(k,2.5k.^(-3),label=L"k^{-3}")
plot!(k,200*k,label=L"k")
ylims!(5e-4,1e2)
vline!([kL],ls=:dash,label=L"k_L")
vline!([kd],ls=:dash,label=L"k_d")
vline!([kξ],ls=:dash,label=L"k_\xi")
xlabel!(L"k\xi")
ylabel!(L"\epsilon_I(k)")