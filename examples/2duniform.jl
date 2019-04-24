using Plots, LaTeXStrings, Pkg, Revise
gr(titlefontsize=12,size=(500,300),transpose=true,colorbar=false)

using FourierGPE

# ==== Units: ========================
# Units of ξ for length, 1/μ for time
# related by μ ≡ ħ²/mξ² = mc^2 where
# c = ħ/mξ is the speed of sound of
# the uniform system.

function showpsi(x,y,ψ)
    p1 = heatmap(x,y,abs2.(ψ),aspectratio=1)
    xlabel!(L"x/\xi");ylabel!(L"y/\xi")
    title!(L"|\psi|^2")
    p2 = heatmap(x,y,angle.(ψ),aspectratio=1)
    xlabel!(L"x/\xi");ylabel!(L"y/\xi")
    title!(L"\textrm{phase}(\psi)")
    p = plot(p1,p2,size=(600,300))
    return p
end

# ==== define user parameters =======
@with_kw mutable struct Params <: UserParams @deftype Float64
    # user parameters:
    κ = 0.1
end
par = Params()

# ==== set simulation parameters ====
L=(200.,200.)
N=(512,512)

# ========= Initialize simulation ======
sim = Sim(L,N,par)
@unpack_Sim sim

# uniform potential
import FourierGPE.V
V(x,y,t) = zero(x*y)

# ========= useful state functions
ψ0(x,y,μ,g) = sqrt(μ/g)*sqrt(max(1.0-V(x,y,0.0)/μ,0.0)+im*0.0)
healing(x,y,μ,g) = 1/sqrt(g*abs2(ψ0(x,y,μ,g)))

x,y = X
#make initial state
ψi = ψ0.(x,y',μ,g)
ϕi = kspace(ψi,sim)
@pack! sim = ϕi
sim

# ====== Evolve in k space ==========
sol = runsim(sim)
# ===================================

# pull out the ground state
ϕg = sol[end]
ψg = xspace(ϕg,sim)
showpsi(x,y,ψg)

# make an initial dipole
using VortexDistributions

ψv = copy(ψg)
d = 10
ξv = healing(0.,0.,μ,g)
dipole = [0.0 d/2 1; 0.0 -d/2 -1]
makeallvortices!(ψv,dipole,x,y,ξv)

# ==== set simulation parameters ====
c = sqrt(μ)
γ = 0.01
tf = 6*L[1]/c

t = LinRange(ti,tf,Nt)
ϕi = kspace(ψv,sim)
reltol = 1e-7 #requires slightly smaller tolerance

@pack! sim = tf,t,γ,ϕi,reltol
# initsim!(sim)

@unpack_Sim sim

# ====== Evolve in k space ==========
solv = runsim(sim)
# ===================================

ψd = xspace(solv[end],sim)
showpsi(x,y,ψd)

# dipole decay
anim = @animate for i=1:Nt
    ψ = xspace(solv[i],sim)
    showpsi(x,y,ψ)
end

gif(anim,"./examples/dipole.gif",fps=30)

# plot energies
function showenergies(ψ,x,y,kx,ky,k2)
    et,ei,ec = energydecomp(ψ,kx,ky',k2)
    p1 = heatmap(x,y,log.(ei),aspectratio=1)
    xlabel!(L"x/\xi");ylabel!(L"y/\xi")
    title!("Incompressible")
    p2 = heatmap(x,y,log.(ec),aspectratio=1)
    xlabel!(L"x/\xi");ylabel!(L"y/\xi")
    title!("Compressible")
    p = plot(p1,p2,size=(600,300))
    return p
end

X,K,dX,dK = makearrays(L,N)
ksq = k2(K)
kx,ky = K
dx,dy = dX

anim = @animate for i=1:Nt
    ψ = xspace(solv[i],sim)
    showenergies(ψ,x,y,K...,k2(K))
end

gif(anim,"./examples/dipoleenergies.gif",fps=30)

# check energy conservation

# nonlinearity
function xenergy(ϕ,sim,t)
    @unpack g,X = sim; x,y = X
    ψ = xspace(ϕ,sim)
    @. ψ *= 0.5*g*abs2(ψ) + V(x,y',t)
    return kspace(ψ,sim)
end

function gpenergy(ϕ,sim,t)
    @unpack μ,γ,espec = sim
    chi = xenergy(ϕ,sim,t)
    H = @. espec*abs2(ϕ) + conj(ϕ)*chi
    return sum(H)*dx*dy |> real
end

H = zero(t)
Ei = zero(t)
Ec = zero(t)
Natoms = zero(t)

for (i,t) in enumerate(t)
    ϕ = solv[i]
    ψ = xspace(ϕ,sim)
    et,ei,ec = energydecomp(ψ,kx,ky',ksq)
    Ei[i] = sum(ei)*dx*dy |> real
    Ec[i] = sum(ec)*dx*dy |> real
    Natoms[i] = sum(abs2.(ψ))*dx*dy
    H[i] = gpenergy(ϕ,sim,t)
end

plot(t,Ei./Natoms,label=L"E_i",legend=:bottomright)
plot!(t,Ec./Natoms,label=L"E_c")
plot!(t,H./Natoms,label=L"H")
xlabel!(L"t")

relerror(x,x0) = (x - x0)/x0 |> abs

plot(t,relerror.(H,H[1]),label = L"relerr(H)",legend=:topleft)
plot!(t,relerror.(Natoms,Natoms[1]),label = L"relerr(N)")


# if we need to save data:
@save "./examples/dipoledecay.jld2" solv.u sim

# another example for saving data as individual files
# using a call back in OrdinaryDiffEq
