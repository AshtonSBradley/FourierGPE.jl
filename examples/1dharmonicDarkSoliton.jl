using Plots, LaTeXStrings 
gr(grid=false,legend=false,titlefontsize=12,size=(500,300),transpose=true,colorbar=false)

using FourierGPE

#--- set simulation parameters ====
L = (40.0,)
N = (512,)
sim = Sim(L,N)
@unpack_Sim sim

μ = 25.0
#--- potential
import FourierGPE.V
V(x,t) = 0.5*x^2

#--- TF state
ψ0(x,μ,g) = sqrt(μ/g)*sqrt(max(1.0-V(x,0.0)/μ,0.0)+im*0.0)

#--- initial state
x = X[1]
ψi = ψ0.(x,μ,g)
ϕi = kspace(ψi,sim)
@pack_Sim! sim

#--- imaginary time to find ground state
sol = runsim(sim)

#--- pull out and check the ground state
ϕg = sol[end]
ψg = xspace(ϕg,sim)
plot(x,one.(x)*μ,ls=:solid,c=:gray,w=1)
plot!(x,V.(x,0.0),c=:red,w=4,alpha=0.4)
plot!(x,g*abs2.(ψg),c=c3,w=2,fill=(0,0.4,c3),size=(600,200))
xlims!(-10,10); ylims!(0,1.3*μ)
title!(L"\textrm{local}\; \mu(x)")
xlabel!(L"x/a_x"); ylabel!(L"\mu(x)/\hbar\omega_x")

#--- imprint dark soliton
ψf = xspace(sol[end],sim)
c = sqrt(μ)
ξ = 1/c
v = 0.5*c
xs = 0.
f = sqrt(1-(v/c)^2)

ψs = ψf.*(f*tanh.(f*(x .-xs)/ξ).+im*v/c);
showpsi(x,ψs)
xlims!(-10,10)

#--- set simulation parameters
γ = 0.0
tf = 8*pi/sqrt(2); t = LinRange(ti,tf,Nt)

# define new simulation using previous parameters
simSoliton = Sim(sim;γ=γ,tf=tf,t=t)
ϕi = kspace(ψs,simSoliton)
@pack! simSoliton = ϕi

#--- evolve in k space
sols = runsim(simSoliton)

#--- plot

ϕf = sols[end-3]
ψf = xspace(ϕf,simSoliton)
showpsi(x,ψf)

anim = @animate for i in 1:length(t)-3
    ψ = xspace(sols[i],simSoliton)
    y = g*abs2.(ψ)
    plot(x,y,c=c3,w=2,fill=(0,0.4,c3),size=(600,300),grid=false)
    xlims!(-10,10); ylims!(0,1.3*μ)
    title!(L"\textrm{local}\; \mu(x)")
    xlabel!(L"x/a_x"); ylabel!(L"\mu(x)/\hbar\omega_x")
end

gif(anim,"./examples/soliton.gif",fps=25)
