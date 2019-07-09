using Test, Plots, LinearAlgebra, Revise, FourierGPE

#2d tests
N = 100

X = xvecs((1,1),(N,N))
K = kvecs((1,1),(N,N))
K2 = k2(K)

ktest = 2*pi
psi = @. exp(im*ktest*X[1]*one.(X[2]'))
# Note: ktest = n*2*pi that are evaluated exactly are the precise values of representation.
# Hence any derivative of a field constructed from a superpoisition of these k's will
# also be exact.

psix = XField(psi,X,K,K2)
psik = KField(fft(psi),X,K,K2)

# flow only in x direction, of correct value
vx,vy = velocity(psix)
@test vx ≈ ktest*one.(vx)
@test vy ≈ zero.(vy)

Vi,Vc = helmholtz(vx,vy,psix)

#orthogonality of Helmholtz decomposition
vidotvc = Vi[1].*Vc[1] .+ Vi[2].*Vc[2]
@test maximum(abs.(vidotvc)) < 1e-10

x,y = X
Vc


# test incompressible
abstract type Vortex end

struct Vortex2 <: Vortex
    xv::Float64
    yv::Float64
    qv::Int64
end

Λ =0.8249
vcore(x) = (Λ*x)^2/(1+(Λ*x)^2)
vdensity(x,y) = vcore(r(x,y))
vphase(x,y) = atan(y,x)
r(x,y) = sqrt(x^2+y^2)
vortexgrid(x,y,x0,y0,q0) = vdensity(x - x0, y - y0)*exp(im*q0*vphase(x-x0,y-y0))

function make_vortex!(psi::XField,vort::Vortex2)
    @unpack xv,yv,qv = vort
    @unpack psiX,X,K,K2=psi; x,y = X
    psiX .*= @. vortexgrid(x,y',xv,yv,qv)
    return XField(psiX,X,K,K2)
end

L = 200
N = 1000

X = xvecs((L,L),(N,N))
K = kvecs((L,L),(N,N))
K2 = k2(K)
psi = one.(X[1].*X[2]') |> complex
psix = XField(psi,X,K,K2)

d = 10.0
vort1 = Vortex2(0.,d/2,1)
vort2 = Vortex2(0.,-d/2,-1)
make_vortex!(psix,vort1)
make_vortex!(psix,vort2)
# x,y=X
# x0 = 0.0; y0 = d/2
# @. psi *= vortexDensityAnsatz(R(x - x0,y' - y0))*exp(im*vortexPhase(x - x0,y' - y0,0.,d/2)-im*vortexPhase(x,y',0.,-d/2))
psi = psix.psiX
rho = abs2.(psi)
heatmap(angle.(psi))
heatmap(rho)

vx,vy = velocity(psix)
jx = rho.*vx
jy = rho.*vy
j = sqrt.(jx.^2 .+ jy.^2)
wx = sqrt.(rho).*vx
wy = sqrt.(rho).*vy
@time Vi,Vc = helmholtz(wx,wy,psix)

ei = 0.5*abs.(Vi[1].^2 .+ Vi[2].^2)
ec = 0.5*abs.(Vc[1].^2 .+ Vc[2].^2)

heatmap(ei)
heatmap(ec)

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
jx = rho.*vx
jy = rho.*vy
j = sqrt.(jx.^2 .+ jy.^2)
wx = sqrt.(rho).*vx
wy = sqrt.(rho).*vy
@time Vi,Vc = helmholtz(wx,wy,psix)

ei = 0.5*abs.(Vi[1].^2 .+ Vi[2].^2)
ec = 0.5*abs.(Vc[1].^2 .+ Vc[2].^2)

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
