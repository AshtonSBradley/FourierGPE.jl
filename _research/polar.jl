using Test, Plots, Parameters, LaTeXStrings
# transform to cartesian => polar coordinates for angular integrals
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2858241/

Nx = 200
xm = 10.
x = LinRange(-xm,xm,Nx); y = x

# create some test data
psi = @. 10*exp(-.1*(x^2+y'^2)) + 15*exp(-(x^2+(y-5)'^2)) # + randn()+im*randn()
heatmap(x,y,abs2.(psi),transpose=true)
xlabel!(L"x");ylabel!(L"y")

abstract type Basis end
struct Oscillator <: Basis
    n::Int64
    ω::Float64
end

Oscillator(n::Int64) = Oscillator(n,1.0)
index(b::Oscillator) = 1:b.n
spectrum(b::Oscillator) = (0:b.n-1)*b.ω
qnumbers(b::Oscillator) = 0:b.n-1

function hermite(x,n)
    @assert n >= 1
    h1 = exp(-x^2 / 2) * (1 / π)^(1 / 4)
    n == 1 && (return h1)
    h2 = sqrt(2) * x * h1
    n == 2 && (return h2)
    h3 = 0.0
    for j ∈ 3:n
        h3 = sqrt(2/(j - 1)) * x * h2 - sqrt((j-2)/(j-1)) * h1
        h1 = h2; h2 = h3
    end
    n > 2 && (return h3)
end
hermite(x,n,ω) = hermite(√ω * x,n) * ω ^ ( 1 / 4 )

function (b::Oscillator)(x::Array{Float64,1})
    @unpack ω,n = b
    T = x * ones(1, n) |> zero
    T[:,1] = @. hermite(x,1,ω)
    n >= 2 && (T[:,2] = @. sqrt(2 * ω) * x * T[:,1])
    for j ∈ index(b)[3:end]
        m = qnumbers(b)[j]
        T[:,j] = @. sqrt(2 / m) * (√ω * x) * T[:,j-1] - sqrt((m-1) / m) * T[:,j-2]
    end
    return T
end

(b::Oscillator)(x::LinRange{Float64}) = b(x |> collect)
(b::Oscillator)(x::Float64,n::Int64) = hermite(x,n,b.ω)

function hermite_polar(x,n,ω)
    T = zeros(size(x)...,n)
    T[:,:,1] = @. hermite(x,1,ω)
    n >= 2 && (T[:,:,2] = @. sqrt(2 * ω) * x * T[:,:,1])
    for j in 3:n
        T[:,:,j] = @. sqrt(2 / (j-1)) * (√ω * x) * T[:,:,j-1] - sqrt((j-2) / (j-1)) * T[:,:,j-2]
    end
    return T
end
hermite_polar(x,n) = hermite_polar(x,n,1.)

# test that we can create the basis, or any one of oscillator modes
N = 30
Nv = 1:N
b1 = Oscillator(N,2.)
@time b1.(x,Nv')
plot(x,hermite.(x,Nv'),legend=false,grid=false)
n = 30
b1 = Oscillator(n)
@time b1(x)
plot(x,b1(x),legend=false,grid=false)

# broadcast over any input parameters
plot(x,hermite.(x,2,[1. .5 .2 .1]),legend=false,grid=false)

# filtering cartesian data
function filterH(psi,x,N)
    H = Oscillator(N)(x)
    T = inv(H'*H)
    F̂ = T*H'*psi*H*T
    psiF = zero(psi)
        for m = 1:N, n = 1:(N+1-m)
            @. psiF += H[:,m]*F̂[m,n]*H[:,n]'
        end
    return psiF
end


# filter data by projecting onto n oscillator modes
n = 30
@time psiF = filterH(psi,x,n)
heatmap(x,y,abs2.(psiF),transpose=true)
xlabel!(L"x");ylabel!(L"y")

# coversion to polar coords
function slowpolar(psi,x,N)
    H = Oscillator(N)(x)
    T = inv(H'*H)
    F̂ = T*H'*psi*H*T

    Nx = length(x)
    θ = LinRange(0,2*pi,2*Nx)
    r = LinRange(0,last(x),Nx/2 |> Int)'

    psiP = zero(θ*r)
    for m = 1:N, n = 1:(N + 1 - m)
        @. psiP += hermite(r*cos(θ),m) * F̂[m,n] * hermite(r*sin(θ),n)
    end
    return psiP
end

function init_polar(x,N)
    θ = LinRange(0,2*pi,2*Nx)
    r = LinRange(0,last(x),Nx/2 |> Int)'
    hx = hermite_polar((@. r*cos(θ)),N)
    hy = hermite_polar((@. r*sin(θ)),N)
    return hx,hy
end

@time hx,hy = init_polar(x,n)

function polar(psi,x,hx,hy)
    N = size(hx)[3]
    H = Oscillator(N)(x)
    T = inv(H'*H)
    F̂ = T*H'*psi*H*T

    Nx = length(x)
    θ = LinRange(0,2*pi,2*Nx)
    r = LinRange(0,last(x),Nx/2 |> Int)'

    psiFP = zero(θ*r)
    for m = 1:N, n = 1:(N + 1 - m)
        Hx = @view hx[:,:,m]
        Hy = @view hy[:,:,n]
        @. psiFP += Hx * F̂[m,n] * Hy
    end
    return psiFP
end

# convert to polar coords
r = LinRange(0,last(x),Nx/2 |> Int)'
θ = LinRange(0,2*pi,2*Nx)

@time psiP = slowpolar(psi,x,30)
@time psifp = polar(psi,x,hx,hy)

# plot in polar coordinates
heatmap(θ,r',abs2.(psiP),transpose=true)
xlabel!(L"\theta");ylabel!(L"r")

heatmap(θ,r',abs2.(psifp),transpose=true)
xlabel!(L"\theta");ylabel!(L"r")

dθ = θ[2]-θ[1]
dr = r[2]-r[1]
dx = x[2]-x[1]
psiPth = sum(psiP,dims=1)*dθ
plot(r',psiPth'.*r')
xlabel!(L"r")

# check total integral preserved (< 1%)
sum(psiPth.*r)*dr
sum(psiF)*dx^2

# ===================================
# ==== Now a more realistic example for trapped BEC
using FourierGPE, VortexDistributions
gr(titlefontsize=12,size=(500,300),transpose=true,colorbar=false)

# ==== Initialize simulation
L = (25.0,25.0)
N = (256,256)
sim = Sim(L,N)
@unpack_Sim sim

# ==== set simulation parameters
μ = 25.0

# Time dependent potential function (here trivial t dep)
import FourierGPE.V
V(x,y,t) = 0.5*(x^2 + y^2)
ψ0(x,y,μ,g) = sqrt(μ/g)*sqrt(max(1.0-V(x,y,0.0)/μ,0.0)+im*0.0)

# ==== make initial state
x,y = X
ψi = ψ0.(x,y',μ,g)
ϕi = kspace(ψi,sim)
@pack_Sim! sim

# ==== evolve
@time sol = runsim(sim)

# speed regression after here?
# ==== ground state
ϕg = sol[end]
ψg = xspace(ϕg,sim)
showpsi(x,y,ψg)

# ==== set simulation parameters
γ = 0.0
t = LinRange(ti,tf,Nt)
ϕi = kspace(ψg,sim)
# reltol = 1e-7
# alg = DP5()

# TF parameters
R(w) = sqrt(2*μ/w^2)
R(1)
rv = 3.
healinglength(x,y,μ,g) = 1/sqrt(g*abs2(ψ0(x,y,μ,g)))
ξ0 = healinglength.(0.,0.,μ,g)
ξ = healinglength(rv,0.,μ,g)

# vcore = Exact(VortexDistributions.ψi,ξ)
vcore = Exact(ξ)
pv = PointVortex(rv,0.,1)
vi = ScalarVortex(vcore,pv)
psi = Torus(copy(ψg),x,y)
vortex!(psi,vi)
showpsi(x,y,psi.ψ)
ψi .= psi.ψ
ϕi = kspace(ψi,sim)



# ===== compare with Fetter JLTP 2010
ξ = 1/sqrt(μ)
Rtf = R(1)
Ωm = 3*log(Rtf/ξ/sqrt(2))/2/Rtf^2
Ωv = Ωm/(1-rv^2/Rtf^2)

# or a precession period
Tv = 2*pi/Ωv

ti = 0.; tf = Tv
t = LinRange(ti,tf,Nt)

@pack_Sim! sim

# ==== evolve
solv = runsim(sim)

# ===================================
ϕf = solv[end]
ψf = xspace(ϕf,sim)
showpsi(x,y,ψf)


anim = @animate for i in eachindex(t)
    ψ = xspace(solv[i],sim)
    showpsi(x,y,ψ)
end

gif(anim,"./examples/vortexprecession.gif",fps=25)
