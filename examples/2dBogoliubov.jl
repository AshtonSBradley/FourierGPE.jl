using Plots, LaTeXStrings, Pkg, Revise, FourierGPE
gr(colorbar=false,size=(600,150),legend=false,grid=false,xticks=true,yticks=true,axis=true)

# ==== set simulation parameters
# an example of saving individual data files for each time
L = (60.0,60.)
N = (512,512)
nfiles = true
filename = "test"
savedir = "data"
path = joinpath(@__DIR__,savedir)
sim = Sim(L,N)

sim = Sim(sim,filename=filename,nfiles=nfiles,path=path)
@unpack_Sim sim

μ = 1.0
g = 0.01
γ = 0.0
tf = 2pi
Nt = 150
ti = 0.0
t = LinRange(ti,tf,Nt)



# ====== Initialize simulation ======


# ===================================
# Bogoliubov state
x,y = X
kx,ky = K

f(k) = 1 + 4/k^2
u(k) = 0.5*(sqrt(f(k))+1)/f(k)^(1/4)
v(k) = 0.5*(sqrt(f(k))-1)/f(k)^(1/4)
lam = 0.01
bog(x,k) = u(k)*exp(im*k*x) - conj(v(k))*exp(-im*k*x)
ψb(x,y,k) = sqrt(μ/g)*(complex(one(x)) + lam*bog(x,k))*one(y)

kb = kx[20]
ψi = ψb.(x,y',kb)
ϕi = kspace(ψi,sim)

# Set time evolution and pack
# reltol = 1e-7
alg = Vern7()

@pack_Sim! sim

# ==== Evolve in k space
@time sol = runsim(sim)
# ===================================

z = abs2.(xspace(sol[end-1],sim)); heatmap(x,y,z)

# ==== 2d movie?
anim = @animate for i in eachindex(t)
    ψ = xspace(sol[i],sim)
    z = abs2.(ψ)
    heatmap(x,y,z,aspectratio=1)
    xlims!(-L[1]/2,L[1]/2)
    ylims!(-L[2]/2,L[2]/2)
end

gif(anim, "./examples/2dBogoliubov.gif", fps = 25)

# ==== slice movie
anim = @animate for i in eachindex(t)
    z = abs2.(xspace(sol[i],sim))[:,256]
    plot(x,z)
    xlims!(-L[1]/2,L[1]/2)
    xlabel!(L"x/\xi")
    ylabel!(L"n(x)")
end

gif(anim, "./examples/2dBogoliubovSlice.gif", fps = 25)


# ==== energy densities

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
    ψ = xspace(sol[i],sim)
    psi = XField(ψ,X,K,K2)
    showenergies(psi)
end

gif(anim,"./examples/2dbogenergies.gif",fps=30)

# energy totals
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
dx = diff(x)[1]; dy = diff(y)[1]

for (i,t) in enumerate(t)
    ϕ = sol[i]
    ψ = xspace(ϕ,sim)
    psi = XField(ψ,X,K,K2)
    et,ei,ec = energydecomp(psi)
    Ei[i] = sum(ei)*dx*dy |> real
    Ec[i] = sum(ec)*dx*dy |> real
    Natoms[i] = sum(abs2.(ψ))*dx*dy
    H[i] = gpenergy(ϕ,sim,t)
end

# ==== compressible energy conservation
plot(t,Ei./Natoms,label=L"E_i",legend=:bottomright)
plot!(t,Ec./Natoms,label=L"E_c")
# plot!(t,H./Natoms,label=L"H")
# xlabel!(L"t")
