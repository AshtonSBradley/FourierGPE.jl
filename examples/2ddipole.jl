#--- load packages
using Test, Plots, LaTeXStrings, Revise, FourierGPE, VortexDistributions
gr(titlefontsize=12,size=(500,300),transpose=true,colorbar=false)

using ColorSchemes
c1 = cgrad(ColorSchemes.linear_blue_5_95_c73_n256.colors)
c2 = cgrad(ColorSchemes.turbo.colors)
import FourierGPE:showpsi

function showpsi(x,y,ψ)
    p1 = heatmap(x,y,abs2.(ψ),aspectratio=1,c=c1,titlefontsize=12,transpose=true,colorbar=false)
    xlims!(x[1],x[end]);ylims!(y[1],y[end])
    xlabel!(L"x");ylabel!(L"y")
    title!(L"|\psi|^2")
    p2 = heatmap(x,y,angle.(ψ),aspectratio=1,c=c2,titlefontsize=12,transpose=true,colorbar=false)
    xlims!(x[1],x[end]);ylims!(y[1],y[end])
    xlabel!(L"x");ylabel!(L"y")
    title!(L"\textrm{phase} (\psi)")
    p = plot(p1,p2,size=(600,300))
    return p
end

function psimovie(sol,sim)
    @unpack X,t = sim; x,y = X
    anim = @animate for i in eachindex(t)
        ψ = xspace(sol[i],sim)
        showpsi(x,y,ψ)
    end
    return anim
end
#---

#--- Set up parameters and ground state
# Units of ξ for length, 1/μ for time
# related by μ ≡ ħ²/mξ² = mc^2 where
# c = ħ/mξ is the speed of sound of
# the uniform system.

# Initialize simulation
L = (200.,200.)
N = (512,512)
sim = Sim(L,N)
@unpack_Sim sim

# set simulation parameters
μ = 1.0
g = 0.01
γ = 0.5
ti = 0.0
tf = .3pi
Nt = 150
t = LinRange(ti,tf,Nt)

# Useful state functions
ψ0(x,y,μ,g) = sqrt(μ/g)*sqrt(max(1.0-V(x,y,0.0)/μ,0.0)+im*0.0)
healinglength(x,y,μ,g) = 1/sqrt(g*abs2(ψ0(x,y,μ,g)))

# Make initial state
x,y = X
kx,ky = K
ψi = ψ0.(x,y',μ,g)
ϕi = kspace(ψi,sim)
@pack_Sim! sim

# evolve
sol = runsim(sim)

# pull out the ground state
ϕg = sol[end]
ψg = xspace(ϕg,sim)
showpsi(x,y,ψg)
#---

#--- periodic dipole phase
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
#---


#--- initial dipole, with periodic phase
ψv = copy(ψg)
d = 14
ξv = healinglength(0.,0.,μ,g)
# can use default ξ = 1
xp,yp = -50.,d/2
xn,yn = -50.,-d/2
pv = PointVortex(xp,yp,1)
nv = PointVortex(xn,yn,-1)

dipole = [nv;pv]
psi = Torus(ψv,x,y) # set field topology for VortexDistributions
vortex!(psi,dipole) # make dipole
ψv = abs.(psi.ψ).*exp.(im*Thetad(x,y,xp,yp,xn,yn))
showpsi(x,y,ψv)
#---

#--- make sure phase is periodic
using Test
testphase = angle.(ψv)
@test testphase[:,1] ≈ testphase[:,end]
@test testphase[1,:] ≈ testphase[end,:]
#---


#--- evolve dipole with weak damping
# Set simulation parameters
c = sqrt(μ)
tf = L[1]/c/2 # a short run for testing evolutoin and detection
γ = 0.03

t = LinRange(ti,tf,Nt)
ϕi = kspace(ψv,sim)
alg = Vern7()

@pack_Sim! sim

# remove boundary artifacts with small γ
solv = runsim(sim)
#---
#--- plot and animate
ψd = xspace(solv[end],sim)
showpsi(x,y,ψd)

anim = psimovie(solv,sim)
gif(anim,"./examples/dipole_damping.gif",fps=25)
#---
#--- Vortex detection test
Nt = length(t)
d = zeros(Nt)

for i=1:Nt
    psi = Torus(xspace(solv[i],sim),x,y)
    vort = findvortices(psi) |> rawData
    d[i] = abs(vort[2,2]-vort[1,2])
end

plot(t[1:Nt],d)
#---
#--- Hamiltonian evolution
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
# ==== simple breakdown
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


# if we need to save data:
# @save "./examples/dipoledecay.jld2" solv.u sim

# another example for saving data as individual files
# using a call back in OrdinaryDiffEq
