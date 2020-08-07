using VortexDistributions, FourierGPE, Plots, LaTeXStrings

## Initialize simulation
# units of healing length, chemical potential
L = (100.,100.)
N = (256,256)
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
# periodic phase correction
# ψd = abs.(psi.ψ).*exp.(im*thetad(x,y,xp,yp,xn,yn))
showpsi(x,y,ψd)

## test new methods
kd = 2*pi/d
kL = 2*pi/L[1]
kξ = 2*pi
kmin = 0.1*kL
kmax = kξ

Np = 400
k = log10range(kmin,kmax,Np)

@time Ek = kespectrum(k,ψd,X,K)
plot(k,Ek,scale=:log10)

## incompressible spectrum
@time Eki = ikespectrum(k,ψd,X,K)
plot!(k,Eki,scale=:log10)

## power-law plot in logspace
n0 = abs2.(ψd[1,1])
E0 = π*n0
p1 = plot(k,Eki/E0,scale=:log10,grid=false,label=L"E_i(k)",legend=:bottomleft)
plot!(k,2k.^(-3),label=L"k^{-3}")
plot!(k,200*k,label=L"k")
ylims!(2e-3,1e2)
# xlims!(0.02,10)
vline!([kL],ls=:dash,label=L"k_L")
vline!([kd],ls=:dash,label=L"k_d")
vline!([kξ],ls=:dash,label=L"k_\xi")
xlabel!(L"k\xi")
ylabel!(L"E_i(k)")

## Analytic form
using SpecialFunctions
vort = findvortices(Torus(ψd,x,y)) |> vortex_array
d = vort[:,2] |> diff
Λ = 0.8249
f(x) = x*(besselk(1,x)*besseli(0,x)-besselk(0,x)*besseli(1,x))
F(x,Λ) = f(x/(2Λ))^2/x
F(x) = F(x,Λ)
Ed(k,d) = 2*F(k)*(1-besselj0(k*d))
plot!(k,Ed.(k,d),label=L"E_i^a(k)")
