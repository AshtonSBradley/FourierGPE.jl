function mixedparsevaltestxk(L,N,j)
X,K,dX,dK = makearrays(L,N)
T = makeT(X,K)

ψ = randn(N...) + im*randn(N...)

n1 = sum(abs2.(ψ))*prod(dX)

ϕ = T.Mxk[j]*ψ

n2 = sum(abs2.(ϕ))*prod(dX)*dK[j]/dX[j]
return n1,n2
end

function mixedparsevaltestkx(L,N,j)
X,K,dX,dK = makearrays(L,N)
T = makeT(X,K)

ϕ = randn(N...) + im*randn(N...)

n1 = sum(abs2.(ϕ))*prod(dK)

ψ = T.Mkx[j]*ϕ

n2 = sum(abs2.(ψ))*prod(dK)*dX[j]/dK[j]
return n1,n2
end

L = (88.)
N = (128)
n1,n2 = mixedparsevaltestxk(L,N,1)
@test n1 ≈ n2

L = (30.,20.)
N = (30,34)
n1,n2 = mixedparsevaltestxk(L,N,1)
@test n1 ≈ n2

n1,n2 = mixedparsevaltestxk(L,N,2)
@test n1 ≈ n2

L = (30.,20.,10)
N = (30,34,24)
n1,n2 = mixedparsevaltestxk(L,N,1)
@test n1 ≈ n2

n1,n2 = mixedparsevaltestxk(L,N,2)
@test n1 ≈ n2

n1,n2 = mixedparsevaltestxk(L,N,3)
@test n1 ≈ n2

L = (30.,20.,10,50.)
N = (30,34,24,22)
n1,n2 = mixedparsevaltestxk(L,N,1)
@test n1 ≈ n2

n1,n2 = mixedparsevaltestxk(L,N,2)
@test n1 ≈ n2

n1,n2 = mixedparsevaltestxk(L,N,3)
@test n1 ≈ n2

n1,n2 = mixedparsevaltestxk(L,N,4)
@test n1 ≈ n2

L = (88.)
N = (128)
n1,n2 = mixedparsevaltestkx(L,N,1)
@test n1 ≈ n2

L = (30.,20.)
N = (30,34)
n1,n2 = mixedparsevaltestkx(L,N,1)
@test n1 ≈ n2

n1,n2 = mixedparsevaltestkx(L,N,2)
@test n1 ≈ n2

L = (30.,20.,10)
N = (30,34,24)
n1,n2 = mixedparsevaltestkx(L,N,1)
@test n1 ≈ n2

n1,n2 = mixedparsevaltestkx(L,N,2)
@test n1 ≈ n2

n1,n2 = mixedparsevaltestkx(L,N,3)
@test n1 ≈ n2

L = (30.,20.,10,50.)
N = (30,34,24,22)
n1,n2 = mixedparsevaltestkx(L,N,1)
@test n1 ≈ n2

n1,n2 = mixedparsevaltestkx(L,N,2)
@test n1 ≈ n2

n1,n2 = mixedparsevaltestkx(L,N,3)
@test n1 ≈ n2

n1,n2 = mixedparsevaltestkx(L,N,4)
@test n1 ≈ n2
