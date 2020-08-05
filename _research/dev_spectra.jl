using Plots, LaTeXStrings, VortexDistributions
gr(titlefontsize=12,size=(500,300),transpose=true,colorbar=false)

using FourierGPE

## Initialize simulation
L = (200.,200.)
N = (256,256)
sim = Sim(L,N)
@unpack_Sim sim

## set simulation parameters
μ = 1.0
g = 0.01
γ = 0.5
ti = 0.0
tf = .3pi
Nt = 150
t = LinRange(ti,tf,Nt)

## useful state functions
ψ0(x,y,μ,g) = sqrt(μ/g)*sqrt(max(1.0-V(x,y,0.0)/μ,0.0)+im*0.0)
healinglength(x,y,μ,g) = 1/sqrt(g*abs2(ψ0(x,y,μ,g)))

## make initial state
x,y = X
kx,ky = K
ψi = ψ0.(x,y',μ,g)

## ground state
ψg = copy(ψi)
showpsi(x,y,ψg)

## initial dipole
ψv = copy(ψg)
d = 30
ξv = healinglength(0.,0.,μ,g)
pv = PointVortex(0.,d/2,1)
nv = PointVortex(0,-d/2,-1)

dipole = [nv;pv]
psi = Torus(ψv,x,y) # set field topology for VortexDistributions
vortex!(psi,dipole) # make dipole
ψv = psi.ψ # pull out the field
showpsi(x,y,ψv)

## set simulation parameters
c = sqrt(μ)
tf = L[1]/c/4
γ = 0.1

t = LinRange(ti,tf,Nt)
ϕi = kspace(ψv,sim)
reltol = 1e-7 #uniform requires slightly smaller tolerance
# alg = Vern7()
simd = Sim(L,N)
@pack_Sim! simd

## remove edge artifacts with small γ
sold = runsim(simd)

## plot
ψd = xspace(sold[end],sim)
showpsi(x,y,ψd)

anim = psimovie(sold,simd)
gif(anim,"./examples/dipole_damping_small.gif",fps=25)


K2 = k2(K)

## energy totals



# all kinetic
heatmap(x,y,log10.(abs2.(fftshift(ϕ))),aspectratio=1)

## single vortex
ψv = copy(ψg)
ξv = healinglength(0.,0.,μ,g)
pv = PointVortex(0.,0.,1)

psi = Torus(ψv,x,y) 
vortex!(psi,pv) 
ψv = psi.ψ 
showpsi(x,y,ψv)

psi = XField(copy(ψv),X,K,K2)
et,ei,ec = ikespectrum(psi,X,K,Kw)
Ei = fftshift(fft(ei))*prod(DX)

heatmap(x,y,log10.(abs2.(fftshift(ϕ))),aspectratio=1)

showenergiesk(psi)
