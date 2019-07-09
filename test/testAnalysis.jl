using Test, Plots, Revise, FourierGPE

N = 100

X = xvecs((1,1),(N,N))
K = kvecs((1,1),(N,N))
K2 = k2(K)

ktest = 2*pi
psi = @. exp(im*ktest*X[1]*one.(X[2]'))
# Note: ktest = n*2*pi that are evaluated exactly are the precise values of representation.
# Hence any derivative of a field constructed from a superpoisition of these k's will
# also be exact.

#2d test

psix = XField(psi,X,K,K2)
psik = KField(fft(psi),X,K,K2)

vx,vy = velocity(psix)
@test vx ≈ ktest*one.(vx)
@test vy ≈ zero.(vy)

Vi,Vc = helmholtz(vx,vy,psix)

vidotvc = Vi[1].*Vc[1] .+ Vi[2].*Vc[2]

maximum(abs2.(vidotvc))

x,y = X
Vc
# test incompressible

using LinearAlgebra, VortexDistributions

L = 200
N = 1000

X = xvecs((L,L),(N,N))
K = kvecs((L,L),(N,N))
K2 = k2(K)
psi = one.(X[1].*X[2]') |> complex
d = 10.0
dipole = [0.0 d/2 1; 0.0 -d/2 -1]
makeallvortices!(psi,dipole,X[1],X[2],1)
rho = abs2.(psi)
heatmap(angle.(psi))
heatmap(rho)

psix = XField(psi,X,K,K2)
psik = KField(fft(psi),X,K,K2)
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
