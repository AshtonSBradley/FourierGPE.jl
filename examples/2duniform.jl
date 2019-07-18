using Test, Plots, LaTeXStrings
gr(titlefontsize=12,size=(500,300),transpose=true,colorbar=false)

using Revise, FourierGPE
# ==== Units: ========================
# Units of ξ for length, 1/μ for time
# related by μ ≡ ħ²/mξ² = mc^2 where
# c = ħ/mξ is the speed of sound of
# the uniform system.

# ==== set simulation parameters
L=(200.,200.)
N=(512,512)
μ = 1.0
# ==== Initialize simulation
sim = Sim(L,N)
@pack! sim = μ
@unpack_Sim sim

# ==== useful state functions
ψ0(x,y,μ,g) = sqrt(μ/g)*sqrt(max(1.0-V(x,y,0.0)/μ,0.0)+im*0.0)
healinglength(x,y,μ,g) = 1/sqrt(g*abs2(ψ0(x,y,μ,g)))


# ==== make initial state
x,y = X
ψi = ψ0.(x,y',μ,g)
ϕi = kspace(ψi,sim)
@pack! sim = ϕi,μ

# ====== Evolve in k space ==========
sol = runsim(sim)

# ==== pull out the ground state
ϕg = sol[end]
ψg = xspace(ϕg,sim)
showpsi(x,y,ψg)

# make an initial dipole
using VortexDistributions

ψv = copy(ψg)
d = 10
ξv = healinglength(0.,0.,μ,g)
@test ξv == 1.0
dipole = [0.0 d/2 1; 0.0 -d/2 -1]
makeallvortices!(ψv,dipole,x,y,ξv)

# ==== set simulation parameters ====
c = sqrt(μ)
γ = 0.01
tf = L[1]/c

t = LinRange(ti,tf,Nt)
ϕi = kspace(ψv,sim)
reltol = 1e-7 #uniform requires slightly smaller tolerance
alg = Vern7()
@pack! sim = tf,t,γ,ϕi,reltol,alg

# ==== Evolve in k space
solv = runsim(sim)

# ==== plot
@unpack_Sim sim

ψd = xspace(solv[end],sim)
showpsi(x,y,ψd)

# dipole decay
anim = @animate for i=1:Nt
    ψ = xspace(solv[i],sim)
    showpsi(x,y,ψ)
end

gif(anim,"./examples/dipole.gif",fps=30)

# plot energies
function showenergies(ψ)
    et,ei,ec = energydecomp(ψ)
    p1 = heatmap(x,y,log10.(ei),aspectratio=1)
    xlabel!(L"x/\xi");ylabel!(L"y/\xi")
    title!("Incompressible")
    p2 = heatmap(x,y,log10.(ec),aspectratio=1)
    xlabel!(L"x/\xi");ylabel!(L"y/\xi")
    title!("Compressible")
    p = plot(p1,p2,size=(600,300))
    return p
end

K2 = k2(K)

anim = @animate for i=1:Nt
    ψ = xspace(solv[i],sim)
    psi = XField(ψ,X,K,K2)
    showenergies(psi)
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
    psi = XField(ψ,X,K,K2)
    et,ei,ec = energydecomp(psi)
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
# @save "./examples/dipoledecay.jld2" solv.u sim

# another example for saving data as individual files
# using a call back in OrdinaryDiffEq
