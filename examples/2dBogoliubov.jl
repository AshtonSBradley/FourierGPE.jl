using Plots, LaTeXStrings
gr(colorbar=false,size=(300,300),legend=false,grid=false)

using FourierGPE

## setup
L = (60.0,60.)
N = (512,512)
sim = Sim(L,N)
@unpack_Sim sim

# individual data files for each time
# nfiles = true
# filename = "test"
# savedir = "data"
# path = joinpath(@__DIR__,savedir)

μ = 1.0
g = 0.01
γ = 0.0
tf = 2pi
Nt = 150
ti = 0.0
t = LinRange(ti,tf,Nt)

## Bogoliubov state
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
# alg = Vern7()

@pack_Sim! sim

## Evolve in k space
@time sol = runsim(sim)

z = abs2.(xspace(sol[end-1],sim)); heatmap(x,y,z)

## movie
anim = @animate for i in eachindex(t)
    ψ = xspace(sol[i],sim)
    n = abs2.(ψ)
    heatmap(x,y,n,aspectratio=1,transpose=true)
    xlims!(-L[1]/2,L[1]/2)
    ylims!(-L[2]/2,L[2]/2)
    xlabel!(L"x/\xi"); ylabel!(L"y/\xi")
end

gif(anim, "./examples/2dBogoliubov.gif", fps = 25)

## slice movie
anim = @animate for i in eachindex(t)
    z = abs2.(xspace(sol[i],sim))[:,256]
    plot(x,z)
    xlims!(-L[1]/2,L[1]/2)
    xlabel!(L"x/\xi")
    ylabel!(L"n(x)")
end

gif(anim, "./examples/2dBogoliubovSlice.gif", fps = 25)


## energy densities

function showenergies(ψ)
    et,ei,ec = energydecomp(ψ)
    p1 = heatmap(x,y,log10.(ei),aspectratio=1,transpose=true)
    xlabel!(L"x/\xi");ylabel!(L"y/\xi")
    title!("Incompressible")
    p2 = heatmap(x,y,log10.(ec),aspectratio=1,transpose=true)
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

gif(anim,"./examples/2dbogenergies.gif",fps = 25)

## energy totals
function xenergy(ϕ,sim,t)
    @unpack g,X = sim; x,y = X
    ψ = xspace(ϕ,sim)
    @. ψ *= 0.5*g*abs2(ψ) + V(x,y',t)
    return kspace(ψ,sim)
end

function gpenergy(ϕ,sim,t)
    @unpack μ,γ,espec,K = sim; kx,ky = K
    dkx,dky = kx[2]-kx[1],ky[2]-ky[1]
    χ = xenergy(ϕ,sim,t)
    H = @. espec*abs2(ϕ) + conj(ϕ)*χ
    return sum(H)*dkx*dky |> real
end

function qpressure(ϕ,sim)
    @unpack g,X,K,espec = sim; x,y = X; kx,ky = K
    dkx,dky = kx[2]-kx[1],ky[2]-ky[1]
    ψ = xspace(ϕ,sim)
    sqrtn = abs.(ψ)
    sqrtnk = kspace(sqrtn,sim)
    Eqp = @. espec*abs2(sqrtnk)
    return sum(Eqp)*dkx*dky |> real
end

H = zero(t)
Ex = zero(t)
Ei = zero(t)
Ec = zero(t)
Eqp = zero(t)
Et = zero(t)
Natoms = zero(t)
dx = diff(x)[1]; dy = diff(y)[1]
dkx = diff(kx)[1]; dky = diff(ky)[1]
for (i,t) in enumerate(t)
    ϕ = sol[i]
    ψ = xspace(ϕ,sim)
    psi = XField(ψ,X,K,K2)
    et,ei,ec = energydecomp(psi)
    Ex[i] = sum(conj(ϕ).*xenergy(ϕ,sim,t))*dkx*dky |> real
    Ei[i] = sum(ei)*dx*dy |> real
    Ec[i] = sum(ec)*dx*dy |> real
    Eqp[i] = qpressure(ϕ,sim)
    Et[i] = Ei[i] + Ec[i] + Eqp[i] + Ex[i]
    Natoms[i] = sum(abs2.(ψ))*dx*dy
    H[i] = gpenergy(ϕ,sim,t)
end

## compressible energy conservation
plot(t,Et./Natoms,label=L"E_t",size=(400,200),legend=true)
plot!(t,Ei./Natoms,label=L"E_i",legend=:bottomright)
plot!(t,Ec./Natoms,label=L"E_c")
plot!(t,H./Natoms,label=L"H")
plot!(t,Eqp./Natoms,label=L"E_{qp}")
plot!(t,Ex./Natoms,label=L"E_{x}")
xlabel!(L"t")
