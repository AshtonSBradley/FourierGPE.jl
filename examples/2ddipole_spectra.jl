using Test, Plots, LaTeXStrings, Revise, FourierGPE, VortexDistributions
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
L = (200.,200.)
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
kx,ky = K
ψi = ψ0.(x,y',μ,g)
ϕi = kspace(ψi,sim)
@pack_Sim! sim

# ====== evolve
sol = runsim(sim)

# ==== pull out the ground state
ϕg = sol[end]
ψg = xspace(ϕg,sim)
showpsi(x,y,ψg)

# ==== initial dipole
ψv = copy(ψg)
d = 30
ξv = healinglength(0.,0.,μ,g)
# can use default ξ = 1
pv = PointVortex(0.,d/2,1)
nv = PointVortex(0,-d/2,-1)

dipole = [nv;pv]
psi = Torus(ψv,x,y) # set field topology for VortexDistributions
vortex!(psi,dipole) # make dipole
ψv = psi.ψ # pull out the field
showpsi(x,y,ψv)

# ==== set simulation parameters
c = sqrt(μ)
tf = L[1]/c/4
γ = 0.1

t = LinRange(ti,tf,Nt)
ϕi = kspace(ψv,sim)
reltol = 1e-7 #uniform requires slightly smaller tolerance
alg = Vern7()
simd = Sim(L,N)
@pack_Sim! simd

# ==== remove boundary artifacts with small γ
sold = runsim(simd)

# ==== plot
ψd = xspace(sold[end],sim)
showpsi(x,y,ψd)

anim = psimovie(sold,sim)
gif(anim,"./examples/dipole_damping_small.gif",fps=25)

# ==== Hamiltonian evolution
γ = 0.0
tf = L[1]/c
t = LinRange(ti,tf,Nt)
ϕi = kspace(ψd,sim)

simh = Sim(L,N)
@pack_Sim! simh

# ==== evolve
solh = runsim(simh)

# ==== plot
ψdh = xspace(solh[end],simh)
showpsi(x,y,ψdh)

anim = psimovie(solh,simh)
gif(anim,"./examples/dipole_hamiltonian_small.gif",fps=25)

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
    ψ = xspace(solh[i],simh)
    psi = XField(ψ,X,K,K2)
    showenergies(psi)
end

gif(anim,"./examples/dipoleenergies_small.gif",fps=30)

# ==== energy totals

function venergy(ϕ,sim,t)
    @unpack g,X = sim; x,y = X
    ψ = xspace(ϕ,sim)
    @. ψ *= V(x,y',t)
    return kspace(ψ,sim)
end

function ienergy(ϕ,sim)
    @unpack g,X = sim; x,y = X
    ψ = xspace(ϕ,sim)
    @. ψ *= 0.5*g*abs2(ψ)
    return kspace(ψ,sim)
end

function kinetic(ϕ,sim)
    @unpack espec,K = sim; kx,ky = K
    dkx,dky = kx[2]-kx[1],ky[2]-ky[1]
    return sum(espec.*abs2.(ϕ))*dkx*dky |> real
end

dx = diff(x)[1]; dy = diff(y)[1]
dkx = diff(kx)[1]; dky = diff(ky)[1]
n0 = μ/g

function energies(sol)
    t = sol.t
    # H = zero(t)
    Ekall = zero(t)
    Ekhy = zero(t)
    Ei = zero(t)
    Ec = zero(t)
    Ev = zero(t)
    Eint = zero(t)
    Eqp = zero(t)
    Et = zero(t)
    Natoms = zero(t)

    for (i,t) in enumerate(t)
    ϕ = sol[i]
    ψ = xspace(ϕ,sim)
    psi = XField(ψ,X,K,K2)
    ek,ei,ec = energydecomp(psi)
    Ekhy[i] = sum(ek)*dx*dy |> real
    Ei[i] = sum(ei)*dx*dy |> real
    Ec[i] = sum(ec)*dx*dy |> real
    Ekall[i] = kinetic(ϕ,sim)

    Ev[i] = sum(@. V(x,y',t)*abs2(ψ))*dx*dy |> real
    Eint[i] = sum(@. g/2*(abs2(ψ) - n0)^2)*dx*dy |> real
    Eqp[i] = Ekall[i] - Ekhy[i]

    Et[i] = Ekall[i] + Ev[i] + Eint[i]
    Natoms[i] = sum(abs2.(ψ))*dx*dy
    end
    return Ei,Ec,Ekhy,Eqp,Ev,Eint,Et,Ekall,Natoms
end

@time Ei,Ec,Ekhy,Eqp,Ev,Eint,Et,Ekall,Natoms = energies(solh)

# ==== plot energies
plot(t,Et./Natoms,label=L"E_t",legend=:left)
plot!(t,Eint./Natoms,label=L"E_{int}")
plot!(t,Ekall./Natoms,label=L"E_{kall}")
ylims!(0,1.1*Et[1]./Natoms[1])
xlabel!(L"t")

# ==== decomposition
plot(t,(Ekall .+ Eint)./Natoms,label=L"E_{kall}+E_{int}",w=5,c=:gray,alpha=0.3)
plot!(t,Et./Natoms,label=L"E_t",legend=:left)
plot!(t,Ei./Natoms,label=L"E_k^i")
plot!(t,Ec./Natoms,label=L"E_k^c")
plot!(t,Eint./Natoms,label=L"E_{int}")
plot!(t,Eqp./Natoms,label=L"E_{qp}")
plot!(t,Ekall./Natoms,label=L"E_{kall}")
xlabel!(L"t")

# energy spectrum in k space
dx
dk = diff(K[1])[1]
DX,DK = dfftall(X,K)

ϕ = solh[end]
ψ = xspace(ϕ,simh)
psi = XField(copy(ψ),X,K,K2)
et,ei,ec = energydecomp(psi)
Ei = fftshift(fft(ei))*prod(DX)

sum(Ei)*dk^2 |> real
sum(ei)*dx^2 |> real

function showenergiesk(ψ)
    et,ei,ec = energydecomp(ψ)
    Et = fftshift(fft(et))*prod(DX)
    Ei = fftshift(fft(ei))*prod(DX)
    Ec = fftshift(fft(ec))*prod(DX)
    p1 = heatmap(x,y,log10.(abs.(Ei)),aspectratio=1)
    xlabel!(L"x/\xi");ylabel!(L"y/\xi")
    title!("Incompressible")
    p2 = heatmap(x,y,log10.(abs.(Ec)),aspectratio=1)
    xlabel!(L"x/\xi");ylabel!(L"y/\xi")
    title!("Compressible")
    p = plot(p1,p2,size=(600,300))
    return p
end

showenergiesk(psi)

# all kinetic
heatmap(x,y,log10.(abs2.(fftshift(ϕ))),aspectratio=1)

# ==== single vortex
ψv = copy(ψg)
ξv = healinglength(0.,0.,μ,g)
pv = PointVortex(0.,0.,1)

psi = Torus(ψv,x,y) # set field topology for VortexDistributions
vortex!(psi,pv) # make dipole
ψv = psi.ψ # pull out the field
showpsi(x,y,ψv)

psi = XField(copy(ψv),X,K,K2)
et,ei,ec = energydecomp(psi)
Ei = fftshift(fft(ei))*prod(DX)

heatmap(x,y,log10.(abs2.(fftshift(ϕ))),aspectratio=1)

showenergiesk(psi)
