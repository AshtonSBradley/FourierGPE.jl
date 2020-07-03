using Revise, FourierGPE
using Test

## Initialize simulation
L = (12.0,12.0)
N = (256,256)
sim = Sim(L,N)
@unpack_Sim sim

# set simulation parameters
μ = 12.0

# Time dependent potential function (here trivial t dep)
import FourierGPE.V
V(x,y,t) = 0.5*(x^2 + y^2)
ψ0(x,y,μ,g) = sqrt(μ/g)*sqrt(max(1.0-V(x,y,0.0)/μ,0.0)+im*0.0)

# make initial state
x,y = X
ψi = ψ0.(x,y',μ,g)
ϕi = kspace(ψi,sim)
@pack_Sim! sim
##

## evolve
@time sol = runsim(sim)

# ground state

ϕg = sol[end]
ψg = xspace(ϕg,sim)
showpsi(x,y,ψg)
##


##
using VortexDistributions
R(w) = sqrt(2*μ/w^2)
R(1)
rv = 2.
healinglength(x,y,μ,g) = 1/sqrt(g*abs2(ψ0(x,y,μ,g)))
ξ0 = healinglength.(0.,0.,μ,g)
ξ = healinglength(rv,0.,μ,g)

pv = PointVortex(rv,0.,1)
vi = ScalarVortex(ξ,pv)


psi = Torus(copy(ψg),x,y)
vortex!(psi,vi)
showpsi(x,y,psi.ψ)

ψi .= psi.ψ
ϕi = kspace(ψi,sim) |> fftshift
kx,ky = K .|> fftshift
x,y = X

## Transform onto polar coordinates in k space
# grids


# make transform to polar grid





function cart2pol(ψi,x,y,kx,ky)
    kmax = maximum(kx)
    nx = length(x)
    nk = nx
    Nk = nk/2 |> Int
    Nθ = 2*nk
    k = LinRange(0,kmax,Nk)
    θ = LinRange(0,2*pi,Nθ+1)[1:Nθ]
    xi = reshape(x,1,1,nx)
    ki = reshape(k,1,1,Nk)
    θi = reshape(θ,1,1,Nθ)
    dθ = diff(θ)[1]
    dk = diff(k)[1]
    dx = diff(x)[1]
    dy = diff(y)[1]

    T = zeros(eltype(ψi),Nk,Nθ,nx)
    S = zeros(eltype(ψi),nx,Nk,Nθ)
    @. T = exp(-im*xi*k*cos(θ'))*dx/sqrt(2*pi)
    @. S = exp(-im*y*k'sin(θi))*dy/sqrt(2*pi)

    ϕp = zeros(eltype(ψi),Nk,Nθ)
    for i in eachindex(x), j in eachindex(y)
        t = @view T[:,:,i]
        s = @view S[j,:,:]
        @. ϕp += t*ψi[i,j]*s
    end
    return ϕp
end

##
@time ϕp = cart2pol(ψi,x,y,kx,ky)
sum(abs2.(ϕp).*k)*dθ*dk
sum(abs2.(ψi))*dx*dy

##
heatmap(log.(eps.() .+abs2.(ϕp)))
Ek = sum(abs2.(ϕp).*k.^3*dθ,dims = 2)
plot(k,log.(Ek))
plot!(k,-k .+5)
vline!([2*pi/ξ0])
