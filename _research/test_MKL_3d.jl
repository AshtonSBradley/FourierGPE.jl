#================================#
#=== 3D expansion simulations ===#
#================================#
using Pkg
Pkg.activate(".")
using FourierGPE, Test


#--- density with Makie
using Makie, AbstractPlotting
function dense(phi)
    ψm = xspace(phi,sim)
    density = abs2.(ψm)
    pmax = maximum(density)
    return density/pmax
end

function densityfilm(sol,Nt;file="3dquenchtrap.gif")
    dir = "media"
    !isdir(dir) && mkdir(dir)
    saveto=joinpath("media",file)

    scene = Scene()
    tindex = Makie.Node(1)
    scene = volume(lift(i -> dense(sol[i]), tindex),
    algorithm = :iso,
    color = (:mediumseagreen,0.15),
    show_axis=false,
    isovalue=3f0(.15))

    R = 70
    eyeat = Vec3f0(R,R,R)
    lookat = Vec3f0(0,0,0)

    record(scene, saveto, 1:Nt) do i
        update_cam!(scene, eyeat, lookat)
        rotate_cam!(scene, 0., -0.1, 0.)
        tindex[] = i
    end
    p = scene[end];
    return
end

#==================================#

#--- 3D expansion
using Plots, Interact, LaTeXStrings, FourierGPE
gr(titlefontsize=12,size=(500,300),transpose=true,colorbar=false)

#--- Scaling GPE
import FourierGPE: nlin!, Lgp!

function k2_scaling(K,λ,t)
    Ks = K./λ(t)
    k2(Ks)
end

#TODO we need to learn a bit about iterators, do, map

# const yt = y'
# const zt = reshape(z,1,1,length(z))

function nlin!(dϕ,ϕ,sim::Sim{3},t)
    @unpack g,X,V0 = sim; x,y,z = X
    yt = y'
    zt = reshape(z,1,1,length(z))
    dϕ .= ϕ
    xspace!(dϕ,sim)
    @. dϕ *= V0 + g*abs2(dϕ)/λx(t)/λy(t)/λz(t) + V(λx(t)*x,λy(t)*yt,λz(t)*zt,t)
    kspace!(dϕ,sim)
end

function Lgp!(dϕ,ϕ,sim::Sim{3},t)
    @unpack γ,μ,espec,K = sim
    nlin!(dϕ,ϕ,sim,t)
    espec = 0.5*k2_scaling(K,λ,t)
    @. dϕ = -im*(1.0 - im*γ)*(dϕ + (espec - μ)*ϕ)
end

function r2_scaling(X,λ,λ̇,t)
    @assert eltype(sqrt.(λ(t).*λ̇(t))) <: Real
    rs = X.*sqrt.(λ(t).*λ̇(t))
    ρs = Iterators.product(rs...)
    map(x->sum(abs2.(x)),ρs)
end

function remove_phase!(ψ,t)
    ψ .*= exp.(-im*r2_scaling(X,λ,λ̇,t)/2)
end

function restore_phase!(ψ,t)
    ψ .*= exp.(im*r2_scaling(X,λ,λ̇,t)/2)
end

init_phase!(ψ) = remove_phase!(ψ,0)

#--- First, we verify Castin-Dum scaling inversion
# using parameters of Mewes et al 1996 (or smaller N)

# initialize default sim
L = (50,50,3)
N = (32,32,32)
sim = Sim(L,N)
@unpack_Sim sim;

# convenient units
mum = 1e-6
nm = 1e-9
amu = 1.661e-27
ħ = 1.055e-34
a0 = 0.05nm

# Rb87 parameters
m = 87amu
as = 5.8nm

# simulation params in units of wr, ar
wr = 18*2pi
wz = 320*2pi
az = sqrt(ħ/m/wz)
ar = sqrt(ħ/m/wr)
gdim = 4*pi*ħ^2*as/m
g = gdim/(ħ*wr)/ar^3


ar = 1.
az = wr/wz |> sqrt
wz /= wr
wr = 1.

γ = 0.1
μ = 50.0
Rr_tf = sqrt(2*μ)
Rz_tf = sqrt(2*μ/wz^2)

# make grids
X,K,dX,dK = makearrays(L,N)
x,y,z = X
yt = y'
zt = reshape(z,1,1,length(z))

# time info
tf = .2/γ
Nt = 100
t = LinRange(0.,tf,Nt)

# potential for initial state (prolate)
import FourierGPE:V
V(x,y,z,t) = 0.5*(x^2 + y^2 + wz^2*z^2)
λx(t) = 1.
λy(t) = 1.
λz(t) = 1.

λẋ(t) = 0.
λẏ(t) = 0.
λż(t) = 0.

λ(t) = (λx(t),λy(t),λz(t))
λ̇(t) = (λẋ(t),λẏ(t),λż(t))

# random initial state
# ψi = randn(N)+im*randn(N)

#Thomas-Fermi initial state
zmid = N[3]/2 |> Int
ψ0(x,y,z,μ,g) = sqrt(μ/g)*sqrt(max(1.0-V(x,y,z,0.0)/μ,0.0)+im*0.0)
ψi = ψ0.(x,yt,zt,μ,g)
showpsi(x,y,ψi[:,:,zmid])

ϕi = kspace(ψi,sim)

# check atom number
Na = sum(abs2.(ψi))*prod(dX)

#save files?
nfiles = true
path = joinpath(@__DIR__,"oblate_ground_state")

#--- Evolve in k space
@pack_Sim! sim
@time sol = runsim(sim)

#--- slice through z=0:
ϕg = sol[end]

# loadfile = joinpath(path,"save200.jld2")
# @load loadfile ψ
# ϕg = ψ
ψg = xspace(ϕg,sim)
showpsi(x,y,ψg[:,:,zmid])

#--- one D plot
Vplot(x,y,z,t) = 0.5*(x^2 + y^2 + wz^2*z^2)
ψtfplot(x,y,z,μ,g)=sqrt(μ/g)*sqrt(max(1.0-Vplot(x,y,z,0.0)/μ,0.0)+im*0.0)
Plots.plot(x,g*abs2.(ψi[:,16,16]),label=L"|\psi(x)|^2")
Plots.plot!(x,g*abs2.(ψtfplot.(x,0,0,μ,g)),label=L"n_{TF}(x)")

#--- file for easy load
function loadpsi(i::Int,sim)
    padi = lpad(string(i),ndigits(length(sim.t)),"0")
    fromfile = joinpath(sim.path,sim.filename*padi*".jld2")
    ϕt = load(fromfile,"ψ")
end


#--- initial state for expansion
# tind = 200
# ϕt = sol[tind]
# ψt = xspace(ϕt,sim)

tind = 100
ϕt = loadpsi(tind,sim)
ψt = xspace(ϕt,sim)

@manipulate for i in 1:zmid
    showpsi(x,y,ψt[:,:,i])
end

#--- save initial state
# using JLD2, FileIO

tofile = "initial_state_oblate.jld2"
save(tofile,"ϕt",ϕt)

#--- Free expansion
fromfile = tofile
load(fromfile,"ϕt")

#save files?
nfiles = true
path = joinpath(@__DIR__,"oblate_expansion")

# potential
V(x,y,z,t) = 0.0

# Scaling GPE: free expansion scaling law (Castin-Dum)
λx(t) = sqrt(1+(wr*t)^2)
λy(t) = sqrt(1+(wr*t)^2)
λz(t) = sqrt(1+(wz*t)^2)

λẋ(t) = wr^2*t/λx(t)
λẏ(t) = wr^2*t/λy(t)
λż(t) = wz^2*t/λz(t)

λ(t) = (λx(t),λy(t),λz(t))
λ̇(t) = (λẋ(t),λẏ(t),λż(t))

# Parameters
γ = 0.0
ti = 0.0
tf = 10
Nt = 200
t = LinRange(ti,tf,Nt)
ϕi .= ϕt

ψi .= xspace(ϕi,sim)
init_phase!(ψi)
ϕi .= kspace(ψi,sim)

simes = Sim(L,N)
@pack_Sim! simes

simes
#--- Evolve
@time soles = runsim(simes)

#--- Plot
@manipulate for i in 1:zmid, tind in 1:length(t)
    ϕf = loadpsi(tind,simes)
    ψf = xspace(ϕf,simes)
    s = t[tind]
    # showpsi(x*λx(s),y*λy(s),ψf[:,:,i])
    # showpsi(x*λy(s),y*λz(s),ψf[i,:,:])
    p1 = Plots.heatmap(y*λy(s),z*λz(s),abs2.(ψf[i,:,:]),c=:inferno,colorbar=true)
    Plots.xlabel!(L"y/a_z")
    Plots.ylabel!(L"z/a_z")
end

#--- aspect ratio
aspect = zero.(t)

x = X[1]
z = X[3]

for i in eachindex(t)
    ϕf = loadpsi(i,simes)
    ψf = xspace(ϕf,simes)
    s = t[i]
    zt = reshape(z,1,1,length(z))
    sx = sum(abs2.(ψf).*x.^2) |> sqrt
    sz = sum(abs2.(ψf).*zt.^2) |> sqrt
    aspect[i] = (λz(s)/λx(s))*sz/sx
end

# aspect_tinf = 2/π/ϵ
# aspect_analytic(t) = (Rr_tf/Rz_tf)*λx(t)/λz(t)

p1 = Plots.plot(t,aspect,size=(500,250),legend=:bottomright,label="sGPE")
# Plots.plot!(t,aspect_analytic.(t),c=:red,label="Castin-Dum")
# Plots.ylims!(1e-1,1e1)
# Plots.plot!(t,one.(t)*aspect_tinf,ls=:dash,label=L"2\omega_\perp/(\pi\omega_z)")
Plots.xlabel!(L"\omega_z t")
Plots.ylabel!(L"R_z(t)/R_\perp(t)")

# <x^2>/<z^2> and Rx^2/Rz^2 are equivalent for Thomas-Fermi approximation

# savefig(p1,"expand_oblate_comparison.pdf")
#--- make movie
# @time p = densityfilm(soles,Nt,file="3dfree_expansion.gif")

#--- plot column density
@manipulate for i in 1:4:lastindex(t)
    ϕf = soles[i]
    ψf = xspace(ϕf,simes)
    ψcd = dropdims(sum(ψf,dims=3),dims=3)
    s = t[i]
    showpsi(x*λ(s),y*λ(s),ψcd/λ(s))
end

#--- compute aspect ratio from standard deviations in r,z
dx = diff(x)[1]; dy = diff(y)[1];dz = diff(z)[1]
Na(ψ) = sum(abs2.(ψ))*dx*dy*dz
Rz(ψ) = sum(abs2.(ψ).*zt.^2)*dx*dy*dz/Na(ψ) |> abs |> sqrt
Rr(ψ) = sum(abs2.(ψ).*(x.^2 .+ yt.^2))*dx*dy*dz/Na(ψ) |> abs |> sqrt
aspect_ratio(ψ) = Rr(ψ)/Rz(ψ)

α0 = aspect_ratio(xspace(soles[1],simes))
α = zero(t)

for j in eachindex(t)
    α[j] = aspect_ratio(xspace(soles[j],simes))
end

plot(t,α,legend=false)
ylims!(0,3.2)
plot!(t,one.(t)*2/α0/pi,ls=:dash)

#--- check against TF radii with fitting
Rz_tf0 = sqrt(2*μ)
Rperp_tf0 = sqrt(2*μ/64)

a0 = Rperp_tf0/Rz_tf0

2/pi/a0
