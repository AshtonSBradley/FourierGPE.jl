using Plots, LaTeXStrings, Pkg, Revise
gr(colorbar=false,size=(600,150),legend=false,grid=false)

using FourierGPE

#--- Initialize simulation
L = (40.0,)
N = (2048,)
sim = Sim(L,N)
@unpack_Sim sim

#--- set simulation parameters
μ = 25.0
γ = 0.0
g = -1.0
n = 5.0
tf = pi/2
Nt = 200
ti = 0.0
t = LinRange(ti,tf,Nt)

#--- soliton wavefunction
x = X[1]
ψs(x) = n*sech(x)
ψi = ψs.(x)
plot(x,abs2.(ψi),lw=10,alpha=0.5)
ϕi = kspace(ψi,sim)

#--- Set all fields
@pack_Sim! sim

#--- Evolve in k space
sol = runsim(sim)

#--- plot
plot(x,abs2.(ψi),lw=3,c=:red)
y = abs2.(xspace(sol[Nt],sim))
plot!(x,y,lw=3,c=c3,fill=(0, 0.4,c3))
xlims!(-5,5)

#--- animate
anim = @animate for i=1:Nt
    ψ = xspace(sol[i],sim)
    y = abs2.(ψ)
    # plot(x,abs2.(ψi),lw=3,c=:red,alpha=0.4)
    plot(x,y,c=c3,fill=(0, 0.4,c3),size=(400,300),xticks=false,yticks=false,axis=false)
    xlims!(-2.2,2.2)
    ylims!(0,220)
end

gif(anim, "./examples/1d_temporal_soliton.gif", fps = 25)
