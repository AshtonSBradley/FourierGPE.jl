# using Test, FourierGPE

# Velocity and Helmholtz tests
N = 100
X = xvecs((1,1),(N,N))
K = kvecs((1,1),(N,N))
K2 = k2(K)

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
