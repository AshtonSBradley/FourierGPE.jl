#transform tests
using Revise

using FourierGPE

import FourierGPE: makeT

# abstract type MixedTransformLibrary end

@with_kw mutable struct TransformsNew{D,N} <: TransformLibrary
    Txk::AbstractFFTs.ScaledPlan{Complex{Float64},FFTW.cFFTWPlan{Complex{Float64},-1,false,D},Float64} = 0.1*plan_fft(crandn_array(D))
    Txk!::AbstractFFTs.ScaledPlan{Complex{Float64},FFTW.cFFTWPlan{Complex{Float64},-1,true,D},Float64} = 0.1*plan_fft!(crandn_array(D))
    Tkx::AbstractFFTs.ScaledPlan{Complex{Float64},FFTW.cFFTWPlan{Complex{Float64},1,false,D},Float64} = 0.1*plan_ifft(crandn_array(D))
    Tkx!::AbstractFFTs.ScaledPlan{Complex{Float64},FFTW.cFFTWPlan{Complex{Float64},1,true,D},Float64} = 0.1*plan_ifft!(crandn_array(D))
    Mxk::Array{AbstractFFTs.ScaledPlan{Complex{Float64},FFTW.cFFTWPlan{Complex{Float64},-1,false,D},Float64},1} = makeTMixed(D)[1]
    Mxk!::Array{AbstractFFTs.ScaledPlan{Complex{Float64},FFTW.cFFTWPlan{Complex{Float64},-1,true,D},Float64},1} = makeTMixed(D)[2]
    Mkx::Array{AbstractFFTs.ScaledPlan{Complex{Float64},FFTW.cFFTWPlan{Complex{Float64},1,false,D},Float64},1} = makeTMixed(D)[3]
    Mkx!::Array{AbstractFFTs.ScaledPlan{Complex{Float64},FFTW.cFFTWPlan{Complex{Float64},1,true,D},Float64},1} = makeTMixed(D)[4]
    psi::ArrayPartition = crandnpartition(D,N)
end

function definetransforms(funcs,args,meas,kwargs)
    trans = []
    for (i,fun) ∈ enumerate(funcs)
        push!(trans, fun(args[i]...,flags=kwargs))
    end
    return meas.*trans
end

function makeT(X,K,flags=FFTW.MEASURE)

    N = [ length(X[i]) for i ∈ eachindex(X) ] |> Tuple
    DX,DK = dfftall(X,K)
    dμx = prod(DX)
    dμk = prod(DK)
    ψtest = ones(N...) |> complex

    trans = (plan_fft,plan_fft!,plan_ifft,plan_ifft!)
    meas = (dμx,dμx,dμk,dμk)
    flags = FFTW.MEASURE
    args = ((ψtest,),(ψtest,),(ψtest,),(ψtest,))
    FFTW.set_num_threads(Sys.CPU_THREADS)

    return definetransforms(trans,args,meas,flags)
end

function makeTMixed(X,K,flags=FFTW.MEASURE)
        # @assert length(X) > 1
        N = [ length(X[i]) for i in eachindex(X) ] |> Tuple
        DX,DK = dfftall(X,K)
        ψtest = ones(N...) |> complex

        args = []
        transxk = []
        transxk! = []
        transkx = []
        transkx! = []
        measxk = []
        measkx = []


        for j = 1:length(N)
            push!(args, (ψtest,j))
            push!(transxk, plan_fft)
            push!(transxk!, plan_fft!)
            push!(transkx, plan_ifft)
            push!(transkx!, plan_ifft!)
            push!(measxk, DX[j])
            push!(measkx, DK[j])

        end

        FFTW.set_num_threads(Sys.CPU_THREADS)
        Mxk = definetransforms(transxk,args,measxk,flags)
        Mxk! = definetransforms(transxk!,args,measxk,flags)
        Mkx = definetransforms(transkx,args,measkx,flags)
        Mkx! = definetransforms(transkx!,args,measkx,flags)

    return Mxk,Mxk!,Mkx,Mkx!
end

function makeTMixed(D::Int)
    L = ones(D) |> Tuple
    N = 4*ones(D) .|> Int |> Tuple
    X,K = makearrays(L,N)
    return makeTMixed(X,K)
end

# test transform build
L = (20.,20.,20.)
N = (64,64,128)
X,K = makearrays(L,N)
t1 = makeT(X,K)
t2 = makeTMixed(X,K)
t2[1]
T = Transforms(t1...)

# seems to be working
T2 = TransformsNew{3,5}()
T3 = TransformsNew{3,2}(t1...,t2...,crandnpartition(3,2))

# to access a scratch field (default is 2x2x...)
T3.psi.x

T3.Mxk


# dims 1, and rest need attention with makeallT
function makeallT(X,K)
    t1 = makeT(X,K)
    t2 = makeTMixed(X,K)
    T = TransformsTest(t1...,t2)
    return T
end

function parsevaltest(L,N)
X,K,dX,dK = makearrays(L,N)
T = makeallT(X,K)

ψ = randn(N...) + im*randn(N...)

n1 = sum(abs2.(ψ))*prod(dX)
ϕ = T.Txk*ψ

n2 = sum(abs2.(ϕ))*prod(dK)
return n1,n2
end

function mixedparsevaltest(L,N,j)
X,K,dX,dK = makearrays(L,N)
T = makeallT(X,K)

ψ = randn(N...) + im*randn(N...)

n1 = sum(abs2.(ψ))*prod(dX)

ϕ = T.Tmixed[1+(j-1)*4 |> Int]*ψ

n2 = sum(abs2.(ϕ))*prod(dX)*dK[j]/dX[j]
return n1,n2
end

T = makeallT(X,K)

L = (100.,200.,300)
N = (50,80,60)
X,K = makearrays(L,N)
t1 = makeT(X,K)



using Test
L = (88.)
N = (128)
n1,n2 = parsevaltest(L,N)
@test n1 ≈ n2

L = (30.,20.)
N = (30,34)
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


# ====================

t4 = makeTMixed(3)
#add
makeTMixed(2)

t1 = makeT(X,K)

test3 = TransformsTest{3}()

M = makeTMixed(X,K)

M[1]

t1 = makeTMixed(X,K)
typeof(t1[1])
typeof(t2.Txk)
t2 = makeT(X,K)

T = Transforms2{3}()

# ok, so we can make the right object, parameterized on system dimension

N = 100
trans = (plan_fft,plan_fft!,plan_ifft,plan_ifft!)
meas = (.3,.1,.1,.1)
ψtest = randn(N,N)+im*randn(N,N)
flags = FFTW.MEASURE

# args = ((ψtest,),(ψtest,),(ψtest,),(ψtest,))
args = (ψtest,)
# ====== simpler approach (?) =====


t1 = deftrans(trans,args,flags)

T = Transforms(t1...)

# ==== give transforms a new life outside Sim =====
using FourierGPE

x = gensym()
eval(:($(x) = Transforms{1}()))
getfield(Main,x)
