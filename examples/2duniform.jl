using Plots, LaTeXStrings, Pkg, Revise
gr(titlefontsize=12,size=(500,300),transpose=true,colorbar=false)

using FourierGPE

# ==== Units: ========================
# Units of ξ for length, 1/μ for time
# related by μ ≡ ħ²/mξ² = mc^2 where
# c = ħ/mξ is the speed of sound of
# the uniform system.

# uniform potential
import FourierGPE.V
V(x,y,t) = zero(x*y)

function showpsi(ψ,x,y)
    p1 = heatmap(x,y,abs2.(ψ),aspectratio=1)
    xlabel!(L"x/\xi");ylabel!(L"y/\xi")
    title!(L"|\psi|^2")
    p2 = heatmap(x,y,angle.(ψ),aspectratio=1)
    xlabel!(L"x/\xi");ylabel!(L"y/\xi")
    title!(L"\textrm{phase}(\psi)")
    p = plot(p1,p2,size=(600,300))
    return p
end

# ==== set simulation parameters ====
sim = Par()
γ = 0.5
Nx = 512; Ny = 512
@pack! sim = γ,Nx,Ny
initsim!(sim)
@unpack_Par sim
# ================================

#make convenient arrays
x,y,kx,ky,dx,dy,dkx,dky = makearrays(Lx,Nx,Ly,Ny)

# useful state functions
ψ0(x,y,μ,g) = sqrt(μ/g) |> complex
healing(x,y,μ,g) = 1/sqrt(g*abs2(ψ0(x,y,μ,g)))

# inital state (uniform)
ψi = ψ0.(x,y',μ,g)

# ====== Evolve in k space ==========
ϕi = kspace(ψi,sim)
sol = runsim(ϕi,sim)
# ===================================

# pull out the ground state
ϕg = sol[end]
ψg = xspace(ϕg,sim)
showpsi(ψg,x,y)

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
tf = 6*Lx/c
t = LinRange(ti,tf,Nt)

@pack! sim = γ,tf,t
initsim!(sim)
@unpack_Par sim
# ==================================

ϕv = kspace(ψv,sim)
solv = runsim(ϕv,sim)
ψd = xspace(solv[200],sim)
showpsi(ψd,x,y)

# dipole decay
anim = @animate for i=1:Nt
    showpsi(xspace(solv[i],sim),x,y)
end

gif(anim,"./dipole.gif",fps=30)

# plot energies
function showenergies(ψ,x,y,kx,ky,k2)
    et,ei,ec = energydecomp(ψ,kx,ky,k2)
    p1 = heatmap(x,y,log.(ei),aspectratio=1)
    xlabel!(L"x/\xi");ylabel!(L"y/\xi")
    title!("Incompressible")
    p2 = heatmap(x,y,log.(ec),aspectratio=1)
    xlabel!(L"x/\xi");ylabel!(L"y/\xi")
    title!("Compressible")
    p = plot(p1,p2,size=(600,300))
    return p
end

anim = @animate for i=1:Nt
    showenergies(xspace(solv[i],sim),x,y,kx,ky,k2)
end

gif(anim,"./examples/dipoleenergies.gif",fps=30)

# check energy conservation

# nonlinearity
function xenergy(ϕ,sim,t)
    @unpack g,x,y = sim
    ψ = xspace(ϕ,sim)
    ψ .*= @. 0.5*g*abs2(ψ) + V(x,y',t)
    return kspace(ψ,sim)
end

function gpenergy(ϕ,sim,t)
    @unpack μ,γ,k2 = sim
    chi = xenergy(ϕ,sim,t)
    Hloc = @. 0.5*k2*abs2.(ϕ) + conj(ϕ)*chi
    return sum(Hloc)*dx*dy |> real
end

H = zero(t)
Ei = zero(t)
Ec = zero(t)
Natoms = zero(t)

for (i,t) in enumerate(t)
    ϕ = solv[i]
    ψ = xspace(ϕ,sim)
    et,ei,ec = energydecomp(ψ,kx,ky,k2)
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
using JLD2, FileIO
@save "./examples/dipoledecay.jld2" solv.u

# another example for saving data as individual files
