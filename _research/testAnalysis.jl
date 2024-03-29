using Test, Plots, LinearAlgebra, Revise, FourierGPE


# Test incompressible without loading VortexDistributions
abstract type Vortex end

struct PointVortex <: Vortex
    xv::Float64
    yv::Float64
    qv::Int64
end

Λ = 0.8249
r(x,y) = sqrt(x^2+y^2)
vcore(r) = (Λ*r)^2/(1+(Λ*r)^2)
vdensity(x,y) = vcore(r(x,y))
vphase(x,y) = atan(y,x)
vortexgrid(x,y,x0,y0,q0) = vdensity(x - x0, y - y0)*exp(im*q0*vphase(x-x0,y-y0))

function make_vortex!(psi::XField,vort::PointVortex)
    @unpack xv,yv,qv = vort
    @unpack psiX,X,K,K2=psi; x,y = X
    psiX .*= @. vortexgrid(x,y',xv,yv,qv)
    return XField(psiX,X,K,K2)
end

L = 500
N = 1000

X = xvecs((L,L),(N,N))
K = kvecs((L,L),(N,N))
K2 = k2(K)
psi = one.(X[1].*X[2]') |> complex
psix = XField(psi,X,K,K2)

d = 10.0
vort1 = PointVortex(0.,d/2,1)
vort2 = PointVortex(0.,-d/2,-1)
@time make_vortex!(psix,vort1)
make_vortex!(psix,vort2)

x,y = X
psi = psix.psiX
rho = abs2.(psi)
heatmap(x,y,angle.(psi))
heatmap(x,y,rho)

vx,vy = velocity(psix)
wx = sqrt.(rho).*vx
wy = sqrt.(rho).*vy
W = wx,wy
@time Wi,Wc = helmholtz(W,psix)

ei = 0.5*abs.(Wi[1].^2 .+ Wi[2].^2)
ec = 0.5*abs.(Wc[1].^2 .+ Wc[2].^2)

x,y = X;dx = diff(x)[1]; dy = diff(y)[1]
sum(ei)*dx*dy
sum(ec)*dx*dy
heatmap(x,y,ei);xlims!(-50,50);ylims!(-50,50)
heatmap(x,y,ec);xlims!(-50,50);ylims!(-50,50)

#TODO test compressible
L = 200
N = 1000

X = xvecs((L,L),(N,N))
K = kvecs((L,L),(N,N))
K2 = k2(K)
psi = one.(X[1].*X[2]') |> complex

ktest = 2*pi
psi .+= @. 0.1*cos(ktest*X[1]*one.(X[2]'))
psix = XField(psi,X,K,K2)

rho = abs2.(psi)
heatmap(angle.(psi))
heatmap(rho)

vx,vy = velocity(psix)
wx = sqrt.(rho).*vx
wy = sqrt.(rho).*vy
@time Wi,Wc = helmholtz(wx,wy,psix)

ei = 0.5*abs.(Wi[1].^2 .+ Wi[2].^2)
ec = 0.5*abs.(Wc[1].^2 .+ Wc[2].^2)

heatmap(ei)
heatmap(ec)

#3d test
X = xvecs((1,2,3),(N,N,N))
K = kvecs((1,2,3),(N,N,N))
K2 = k2(K)
psi = randn(N,N,N) + im*randn(N,N,N)

psix = XField(psi,X,K,K2)
psik = KField(fft(psi),X,K,K2)

@time v = velocity(psix)
@time Vi,Vc = helmholtz(v,psix)
