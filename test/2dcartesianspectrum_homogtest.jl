using Test, SpecialFunctions, VortexDistributions
using LazyArrays, FillArrays
using Revise, FourierGPE


#--- Initialize simulation
# units of healing length, chemical potential
L = (100.,100.)
N = (512,512)
sim = Sim(L,N)
@unpack_Sim sim
#--- set simulation parameters
μ = 1.0
g = 0.01
γ = 0.5
ti = 0.0
tf = .3pi
Nt = 150
t = LinRange(ti,tf,Nt)
#--- Useful state functions
ψ0(x,y,μ,g) = sqrt(μ/g)*sqrt(max(1.0-V(x,y,0.0)/μ,0.0)+im*0.0)
healinglength(x,y,μ,g) = 1/sqrt(g*abs2(ψ0(x,y,μ,g)))
#--- Make initial state
x,y = X
kx,ky = K
ψi = ψ0.(x,y',μ,g)
ϕi = kspace(ψi,sim)
@pack_Sim! sim

# evolve
sol = runsim(sim)

#--- pull out the ground state
ϕg = sol[end]
ψg = xspace(ϕg,sim)
showpsi(x,y,ψg)

#--- periodic dipole phase
# Billam et al, PRL 112, 145301 (2014), Supplemental
# Methods
H(x) = x > 0. ? 1.0 : 0.0
shift(x,xi) = x - xi
tans(x,xk) = tan((shift(x,xk) - π)*0.5)
tanhs(x,xk,j) = tanh((shift(x,xk) + 2*π*j)*0.5)

function kernel(x,y,xp,yp,xn,yn,j)
    return atan(tanhs(y,yn,j)*tans(x,xn)) -
    atan(tanhs(y,yp,j)*tans(x,xp))
end

# Dimensionless form
function θd(x,y,dip)
    vp,vn = dip
    xp,yp,_ = rawData(vp)
    xn,yn,_ = rawData(vn)
    s = 0.0
    for j = -5:5
        s += kernel(x,y,xp,yp,xn,yn,j)
    end
    return s + π*(H(shift(x,xp)) - H(shift(x,xn))) - y*(xp - xn)/(2*π)
end
# arbitrary domains and dipole sizes:
function thetad(x,y,xp,yp,xn,yn)
    s = 0.0
    for j = -5:5
        s += kernel(x,y,xp,yp,xn,yn,j)
    end
    s += π*(H(shift(x,xp)) - H(shift(x,xn))) - y*(xp - xn)/(2*π)
    return s - x*H(abs(yp - yn) - π) + y*H(abs(xp - xn) - π)
end

function Thetad(x,y,xp,yp,xn,yn)
    Lx = x[end]-x[1]
    Ly = y[end]-y[1]
    return @. angle(exp(im*thetad.(x*2*pi/Lx,y'*2*pi/Ly,xp*2*pi/Lx,yp*2*pi/Ly,xn*2*pi/Lx,yn*2*pi/Ly)))
end
#--- initial dipole, with periodic phase
ψv = copy(ψg)
d = 14
ξv = healinglength(0.,0.,μ,g)
# can use default ξ = 1
xp,yp = 0.,d/2
xn,yn = 0.,-d/2
pv = PointVortex(xp,yp,1)
nv = PointVortex(xn,yn,-1)

dipole = [nv;pv]
psi = Torus(ψv,x,y) # set field topology for VortexDistributions
vortex!(psi,dipole) # make dipole
ψv = abs.(psi.ψ).*exp.(im*Thetad(x,y,xp,yp,xn,yn))
showpsi(x,y,ψv)

#--- evolve dipole with weak damping
# Set simulation parameters
c = sqrt(μ)
tf = L[1]/c # a short run for testing evolutoin and detection
γ = 0.03

t = LinRange(ti,tf,Nt)
ϕi = kspace(ψv,sim)
alg = Vern7()

@pack_Sim! sim

# remove boundary artifacts with small γ
@time solv = runsim(sim)

#--- pull out final state
ψd = xspace(solv[end],sim)
showpsi(x,y,ψd)

#--- construct polar spectrum
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

# measures for unpadded fields
dx,dy = diff(x)[1],diff(y)[1]
dkx,dky = diff(kx)[1],diff(ky)[1]
DX,DK = dfftall(X,K)
ϕd = fft(ψd)*prod(DX)
@test sum(abs2.(ψd))*dx*dy ≈ sum(abs2.(ϕd))*dkx*dky


ψc = zeropad(ψd)
ϕc = fft(ψc)
A = ifft(abs2.(ϕc)) |> fftshift

kL = 0.1*2*pi/L[1]
kξ =   2*pi
kmin = kL

kmax = kξ
Np = 200
kp = log10range(kmin,kmax,Np)

Ek = kespectrum(kp,ψi,x,y)
plot(kp,Ek,scale=:log10)

# heatmap(log.(abs2.(A |> fftshift) .+ eps.()))

#--- test incompressible spectrum
k2field = k2(K)
psi = XField(ψd,X,K,k2field)
vx,vy = velocity(psi)
rho = abs2.(ψd)
wx = @. sqrt(rho)*vx; wy = @. sqrt(rho)*vy
Wi, Wc = helmholtz(wx,wy,psi)
wix,wiy = Wi

function ikespectrum(kp,wconv,x,y)
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

#--- convolutions
# measures for zero padded fields
DX,DK = dfftall(X,K)
# no change to measures!

# transforms
Wix = zeropad(wix)
Wiy = zeropad(wiy)
Wixk = fft(Wix)*prod(DX)
Wiyk = fft(Wiy)*prod(DX)
@test sum(abs2.(Wix))*dx*dy ≈ sum(abs2.(Wixk))*(dkx/2)*(dky/2)

# convolutions
Cix = ifft(abs2.(Wixk))*prod(DK) |> fftshift
Ciy = ifft(abs2.(Wiyk))*prod(DK) |> fftshift
Ci = Cix .+ Ciy

kL = 2*pi/L[1]
kξ = 2*pi
kd = 2*pi/d

kmin = 0.1*kL
kmax = 1.3kξ
Np = 400
kp = log10range(kmin,kmax,Np)

Eki = ikespectrum(kp,Ci,x,y)

#--- power law plot in logspace
n0 = abs2.(ψd[1,1])
E0 = π*n0
p1 = plot(kp,Eki/E0,scale=:log10,grid=false,label=L"E_i(k)",legend=:bottomleft)
plot!(kp,2kp.^(-3),label=L"k^{-3}")
plot!(kp,200*kp,label=L"k")
ylims!(5e-4,1e2)
# xlims!(0.02,10)
vline!([kR],ls=:dash,label=L"k_L")
vline!([kd],ls=:dash,label=L"k_d")
vline!([kξ],ls=:dash,label=L"k_\xi")
xlabel!(L"k\xi")
ylabel!(L"E_i(k)")

# test against analytic form to get norm of convolution right
vort = findvortices(Torus(ψd,x,y)) |> rawData
d = vort[:,2] |> diff
Λ = 0.8249
f(x) = x*(besselk(1,x)*besseli(0,x)-besselk(0,x)*besseli(1,x))
F(x,Λ) = f(x/(2Λ))^2/x
F(x) = F(x,Λ)
Ed(k,d) = 2*F(k)*(1-besselj0(k*d))

plot!(kp,Ed.(kp,d),label=L"{\cal E}_i^a(k)")
