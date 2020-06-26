using Plots, LaTeXStrings 
gr(colorbar=false,size=(600,150),legend=false,grid=false,xticks=false,yticks=false,axis=false)

using FourierGPE


## Initialize simulation
L = (60.0,)
N = (512,)
sim = Sim(L,N)
@unpack_Sim sim;

## set simulation parameters
μ = 25.0
g = -0.01
γ = 0.0
Ns = 200
ξs = 2/abs(g)/Ns
us = 20
tf = π |> Float64
Nt = 150
ti = 0.0
t = LinRange(ti,tf,Nt);

## soliton wavefunction
x = X[1]
ψs(x) = sqrt(Ns/2ξs)*sech(x/ξs)*exp(im*us*x)
ψi = ψs.(x)
ϕi = kspace(ψi,sim)

## Set all fields
@pack_Sim! sim

## Evolve in k space
sol = runsim(sim)

## plot
y = abs2.(xspace(sol[Nt-76],sim))
plot(x,y,fill=(0, 0.2))

anim = @animate for i=1:Nt-8
    ψ = xspace(sol[i],sim)
    y = abs2.(ψ)
    plot(x,y,c=c3,fill=(0, 0.4,c3))
end

gif(anim, "./examples/brightsoliton.gif", fps = 25)
