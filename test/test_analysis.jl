## analysis
using Test, FourierGPE

# Velocity and Helmholtz tests
n = 100
L = (1,1)
N = (n,n)
X = xvecs(L,N)
K = kvecs(L,N)
K2 = k2(K);

##
ktest = K[1][2]
psi = @. exp(im*ktest*X[1]*one.(X[2]'))

psix = XField(psi,X,K,K2)

# flow only in x direction, of correct value
vx,vy = velocity(psix)
@test vx ≈ ktest*one.(vx)
@test vy ≈ zero.(vy)

# Decomposition
Vi,Vc = helmholtz(vx,vy,psix)

# Orthogonality of Helmholtz decomposition
vidotvc = Vi[1].*Vc[1] .+ Vi[2].*Vc[2]
@test maximum(abs.(vidotvc)) < 1e-10

# Projective
@test Vi[1] .+ Vc[1] ≈ vx
@test Vi[2] .+ Vc[2] ≈ vy

## autocorrelation
x,y = X; kx,ky = K
dx = x[2]-x[1];dk = kx[2]-kx[1]
 
Natoms = sum(abs2.(psi))*dx^2
AC = autocorrelate(psi,X,K)
Na,Nim = AC[101,101] .|> (real,imag) 

@test Natoms ≈ Na
@test Nim ≈ 0.0 


