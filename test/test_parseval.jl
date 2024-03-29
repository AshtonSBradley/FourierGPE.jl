#transform tests

# test Parseval's theorem for wavefunctions
function parsevaltest(L,N)
    X,K,dX,dK = makearrays(L,N)
    T = makeT(X,K)

    ψ = randn(N...) + im*randn(N...)

    n1 = sum(abs2.(ψ))*prod(dX)
    ϕ = T.Txk*ψ

    n2 = sum(abs2.(ϕ))*prod(dK)
    return n1,n2
end

L = (88.,)
N = (128,)
# X,K,dX,dK = makearrays(L,N)
# T = makeT(X,K) # TODO: norm problem in makeT for 2,3 dims

n1,n2 = parsevaltest(L,N)
@test n1 ≈ n2

## debug offending test
L = (30.,20.)
N = (30,34)

# X,K,dX,dK = makearrays(L,N)
# T = makeT(X,K) # TODO: norm problem in makeT for 2,3 dims

# ψ = randn(N...) + im*randn(N...)

# n1 = sum(abs2.(ψ))*prod(dX)
# ϕ = T.Txk*ψ

# n2 = sum(abs2.(ϕ))*prod(dK)

n1,n2 = parsevaltest(L,N)
@test n1 ≈ n2

L = (30.,20.,10)
N = (30,34,24)
n1,n2 = parsevaltest(L,N)
@test n1 ≈ n2

L = (30.,20.,10,50.)
N = (30,34,24,22)
n1,n2 = parsevaltest(L,N)
@test n1 ≈ n2
