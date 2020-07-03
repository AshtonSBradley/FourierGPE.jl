# Test MKL plans
using Pkg
Pkg.activate(".")
using FourierGPE, Test

#transform tests
# using Test


## test Parseval's theorem for wavefunctions
function parsevaltest(L,N)
X,K,dX,dK = makearrays(L,N)
T = makeT(X,K)

ψ = randn(N...) + im*randn(N...)

n1 = sum(abs2.(ψ))*prod(dX)
ϕ = T.Txk*ψ

n2 = sum(abs2.(ϕ))*prod(dK)
return n1,n2
end
##

L = (88.)
N = (128)
n1,n2 = parsevaltest(L,N)
@test n1 ≈ n2

L = (30.,20.)
N = (30,34)
n1,n2 = parsevaltest(L,N)
@test n1 ≈ n2

## Fails:
L = (30.,20.,10)
N = (30,34,24)
X,K,dX,dK = makearrays(L,N)
# T = makeT(X,K)
FFTW.set_num_threads(Sys.CPU_THREADS)
flags=FFTW.MEASURE
dim = length(X)
N = length.(X)
DX,DK = dfftall(X,K)
dμx = prod(DX)
dμk = prod(DK)
ψtest = ones(N...) |> complex

trans = (plan_fft,plan_fft!,plan_ifft,plan_ifft!)
meas = (dμx,dμx,dμk,dμk)
# flags = FFTW.MEASURE
args = ((ψtest,),(ψtest,),(ψtest,),(ψtest,))

# flags != FFTW.MEASURE && @info "Planning FFTs ..."
Txk,Txk!,Tkx,Tkx! = definetransforms(trans,args,meas,flags)

#Fails:
Mxk,Mxk!,Mkx,Mkx! = makeTMixed(X,K,flags=flags)

# @info "...Plans created."
# return Transforms{dim,j}(Txk,Txk!,Tkx,Tkx!,Mxk,Mxk!,Mkx,Mkx!,crandnpartition(dim,j))


# n1,n2 = parsevaltest(L,N)
# @test n1 ≈ n2

##
L = (30.,20.,10,50.)
N = (30,34,24,22)
n1,n2 = parsevaltest(L,N)
@test n1 ≈ n2
