using FourierGPE
gr(titlefontsize=12,size=(500,300),transpose=true,colorbar=false)

# ==== Units
# this example works in oscillator units

# ==== Initialize simulation

L = (40.0,40.0)
N = (256,256)
sim = Sim(L,N)
@unpack_Sim sim

# ==== simulation parameters
μ = 15.0


# ==== potential
import FourierGPE.V
V(x,y,t) = 0.5*(x^2 + y^2)

# ==== TF state
ψtf(x,y,μ,g) = sqrt(μ/g)*sqrt(max(1.0-V(x,y,0.0)/μ,0.0)+im*0.0)

x,y = X
# ==== make initial state
ψi = ψtf.(x,y',μ,g)
ϕi = kspace(ψi,sim)

@pack_Sim! sim

# ==== evolve
sol = runsim(sim)

# ==== pull out the ground state
ϕg = sol[end]
ψg = xspace(ϕg,sim)
showpsi(x,y,ψg)

R(w) = sqrt(2*μ/w^2)
R(1)
plot(x, abs2.(ψtf.(x,0.,μ,g)))
plot!(x, abs2.(ψg[:,128]))

# ==== free expansion
V(x,y,t) = 0.0

# ==== sim parameters
γ = 0.0
tf = 10.0
t = LinRange(ti,tf,Nt)
ϕi = ϕg

sime = Sim(L,N)
@pack! sime = tf,t,γ,ϕi

# ==== evolve
sole = runsim(sime)

ϕf = sole[end]
ψf = xspace(ϕf,sime)
showpsi(x,y,ψf)

anim = @animate for i=1:Nt-6
    ψ = xspace(sole[i],sime)
    showpsi(x,y,ψ)
end

gif(anim,"./examples/expand.gif",fps=30)
