using Plots, LaTeXStrings
gr(size=(600,300),colorbar=false,legend=false,grid=false)

using FourierGPE

## initialize
L = (60.0,)
N = (512,)
sim = Sim(L,N)
@unpack_Sim sim

## simulation parameters
μ = 1.0
g = 0.01
γ = 0.0
tf = 20.0pi
Nt = 150
ti = 0.0
t = LinRange(ti,tf,Nt)

## Bogoliubov state
x = X[1]
k = K[1]

# Define Bog state
f(k) = 1 + 4/k^2
u(k) = 0.5*(sqrt(f(k))+1)/f(k)^(1/4)
v(k) = 0.5*(sqrt(f(k))-1)/f(k)^(1/4)
lam = 0.01
bog(x,k) = u(k)*exp(im*k*x) - conj(v(k))*exp(-im*k*x)
ψb(x,k) = sqrt(μ/g)*(complex(1) + lam*bog(x,k))

## choose an allowed k for periodic system
kb = k[20]
ψi = ψb.(x,kb)
ϕi = kspace(ψi,sim)

# Set time evolution and pack
# reltol = 1e-7
# alg = Vern7()

@pack_Sim! sim

## Evolve in k space
@time sol = runsim(sim)
y = abs2.(xspace(sol[end-1],sim)); plot(x,y)

## density movie
anim = @animate for i=1:Nt
    ψ = xspace(sol[i],sim)
    y = g*abs2.(ψ)
    plot(x,y,c=c3,fill=(0,0.5,c3))
    plot!(x,one.(x)*μ,c=:black,alpha=0.3)
    ylims!(0,1.1)
    xlims!(first(x),last(x))
    xlabel!(L"x/\xi"); ylabel!(L"gn(x)/\mu")
end

gif(anim, "./examples/bogoliubov_density.gif", fps = 20)
