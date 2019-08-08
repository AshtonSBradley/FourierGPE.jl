using Test, Plots, LaTeXStrings, Revise, FourierGPE
gr(titlefontsize=12,size=(500,300),transpose=true,colorbar=false)

function psimovie(sol,sim)
    @unpack X,t = sim; x,y = X
    anim = @animate for i in eachindex(t)
        ψ = xspace(sol[i],sim)
        showpsi(x,y,ψ)
    end
    return anim
end

# ==== Units
# Units of ξ for length, 1/μ for time
# related by μ ≡ ħ²/mξ² = mc^2 where
# c = ħ/mξ is the speed of sound of
# the uniform system.



# ==== Initialize simulation
L = (400.,400.)
N = (256,256)
sim = Sim(L,N)
@unpack_Sim sim

# ==== set simulation parameters
μ = 1.0
g = 0.01
γ = 0.5
ti = 0.0
tf = .3pi
Nt = 150
t = LinRange(ti,tf,Nt)

# ==== useful state functions
ψ0(x,y,μ,g) = sqrt(μ/g)*sqrt(max(1.0-V(x,y,0.0)/μ,0.0)+im*0.0)
healinglength(x,y,μ,g) = 1/sqrt(g*abs2(ψ0(x,y,μ,g)))

# ==== make initial state
x,y = X
ψi = ψ0.(x,y',μ,g)
ϕi = kspace(ψi,sim)
@pack! sim = ϕi,μ

# ====== evolve
sol = runsim(sim)

# ==== pull out the ground state
ϕg = sol[end]
ψg = xspace(ϕg,sim)
showpsi(x,y,ψg)

# ==== initial dipole
using VortexDistributions

ψv = copy(ψg)
d = 30
ξv = healinglength(0.,0.,μ,g)
# can use default healing length = 1
pv = PointVortex(0.,d/2,1)
nv = PointVortex(0,-d/2,-1)

dipole = [nv;pv]
psi = Torus(ψv,x,y)
vortex!(psi,dipole)
ψv = psi.ψ

# ==== set simulation parameters
c = sqrt(μ)
tf = L[1]/c/2
γ = 0.1

t = LinRange(ti,tf,Nt)
ϕi = kspace(ψv,sim)
reltol = 1e-7 #uniform requires slightly smaller tolerance
alg = Vern7()

@pack_Sim! sim

# ==== evolve
solv = runsim(sim)

# ==== plot
ψd = xspace(solv[end],sim)
showpsi(x,y,ψd)

anim = psimovie(solv,sim)
gif(anim,"./examples/dipole_damping.gif",fps=25)

# ==== Hamiltonian evolution
γ = 0.0
tf = L[1]/c
t = LinRange(ti,tf,Nt)
ϕi = kspace(ψd,sim)
@pack_Sim! sim

# ==== evolve
solh = runsim(sim)

# ==== plot
ψdh = xspace(solh[end],sim)
showpsi(x,y,ψdh)

anim = psimovie(solh,sim)
gif(anim,"./examples/dipole_hamiltonian.gif",fps=25)

# energy densities
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

anim = @animate for i in eachindex(t)
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
    @unpack μ,γ,K,espec = sim; kx,ky=K
    dkx,dky = kx[2]-kx[1],ky[2]-ky[1]
    chi = xenergy(ϕ,sim,t)
    H = @. espec*abs2(ϕ) + conj(ϕ)*chi
    return sum(H)*dkx*dky |> real
end

H = zero(t)
Ei = zero(t)
Ec = zero(t)
Natoms = zero(t)
dx,dy = diff(x)[1],diff(y)[1]

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
xlabel!(L"t")
plot!(t,H./Natoms,label=L"H")


relerror(x,x0) = (x - x0)/x0 |> abs

plot(t,relerror.(H,H[1]),label = L"relerr(H)",legend=:topleft)
plot!(t,relerror.(Natoms,Natoms[1]),label = L"relerr(N)")


# if we need to save data:
# @save "./examples/dipoledecay.jld2" solv.u sim

# another example for saving data as individual files
# using a call back in OrdinaryDiffEq
