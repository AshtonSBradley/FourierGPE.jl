#transform tests
using Revise

using FourierGPE, RecursiveArrayTools

N = 100
trans = (plan_fft,plan_fft!,plan_ifft,plan_ifft!)
meas = (.3,.1,.1,.1)
ψtest = randn(N)+im*randn(N)
flags = FFTW.MEASURE

# args = ((ψtest,),(ψtest,),(ψtest,),(ψtest,))
args = (ψtest,)
# ====== simpler approach (?) =====

function deftrans(funcs,args,kwargs)
    trans = []
    for fun ∈ funcs
        push!(trans, fun(args...,flags=kwargs))
    end
    return meas.*trans
end

t1 = deftrans(trans,args,flags)


# ====== simpler approach =====

T = Transforms(t1...)

# test vecgtor of array
using FastGaussQuadrature, Test

x,w = gausslaguerre(70)
T2 = VectorOfArray([w*w'])
push!(T2,x*x')


recs = [[1, 2, 3], [4, 5, 6], [7, 8, 9]]
testva = VectorOfArray(recs)

for (i, elem) in enumerate(testva)
    @test elem == testva[i]
end

push!(testva, [10, 11, 12])

testva[1]


recs = [[1 2 3;4 5 6], [4 5 6; 7 8 9], [7 8 9; 10 11 12]]
testva = VectorOfArray(recs)

for (i, elem) in enumerate(testva)
    @test elem == testva[i]
end
testva[1]

append!(testva, [10 11 12; 13 14 15] )



# ==== add mixed transforms =====
abstract type MixedTransformLibrary end

@with_kw mutable struct Transforms2{D} <: TransformLibrary
    Txk::AbstractFFTs.ScaledPlan{Complex{Float64},FFTW.cFFTWPlan{Complex{Float64},-1,false,D},Float64} = 0.1*plan_fft(crandn_array(D))
    Txk!::AbstractFFTs.ScaledPlan{Complex{Float64},FFTW.cFFTWPlan{Complex{Float64},-1,true,D},Float64} = 0.1*plan_fft!(crandn_array(D))
    Tkx::AbstractFFTs.ScaledPlan{Complex{Float64},FFTW.cFFTWPlan{Complex{Float64},1,false,D},Float64} = 0.1*plan_ifft(crandn_array(D))
    Tkx!::AbstractFFTs.ScaledPlan{Complex{Float64},FFTW.cFFTWPlan{Complex{Float64},1,true,D},Float64} = 0.1*plan_ifft!(crandn_array(D))
    Tmixed::MixedTransformLibrary = makeTMixed()
end


using FFTW

N = 1024
psitest = randn(N,N) + im*randn(N,N)
P = plan_fft(psitest,2)

phitest = P*psitest

import FourierGPE: makeT

function makeT(X,K,flags=FFTW.MEASURE)

    N = [ length(X[i]) for i ∈ eachindex(X) ] |> Tuple
    DX,DK = dfftall(X,K)
    dμx = prod(DX)
    dμk = prod(DK)
    ψtest = ones(N...) |> complex
    #
    # # plan transforms
    # FFTW.set_num_threads(Sys.CPU_THREADS)
    # Txk = dμx*plan_fft(ψtest,flags=flags)
    # Txk! = dμx*plan_fft!(ψtest,flags=flags)
    # Tkx  = dμk*plan_ifft(ψtest,flags=flags)
    # Tkx!  = dμk*plan_ifft!(ψtest,flags=flags)
    # T = Transforms(Txk,Txk!,Tkx,Tkx!)

    trans = (plan_fft,plan_fft!,plan_ifft,plan_ifft!)
    meas = (dμx,dμx,dμk,dμk)
    flags = FFTW.MEASURE
    FFTW.set_num_threads(Sys.CPU_THREADS)
    args = ((ψtest,),(ψtest,),(ψtest,),(ψtest,))

    t1 = deftrans(trans,args,meas,flags)

    T = Transforms(t1...)

    return T
end

function definetransforms(funcs,args,meas,kwargs)
    trans = []
    for (i,fun) ∈ enumerate(funcs)
        push!(trans, fun(args[i]...,flags=kwargs))
    end
    return meas.*trans
end

L = (100.,200.,300)
N = (50,80,60)
X,K = makearrays(L,N)
T = makeT(X,K)

function parsevaltest(L,N)
X,K,dX,dK = makearrays(L,N)
T = makeT(X,K)

ψ = randn(N...) + im*randn(N...)

n1 = sum(abs2.(ψ))*prod(dX)
ϕ = T.Txk*ψ

n2 = sum(abs2.(ϕ))*prod(dK)
return n1,n2
end

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


function makeTMixed(X,K,flags=FFTW.MEASURE)
        @assert length(X) > 1
        N = [ length(X[i]) for i ∈ eachindex(X) ] |> Tuple
        DX,DK = dfftall(X,K)
        ψtest = ones(N...) |> complex

        flags = FFTW.MEASURE
        FFTW.set_num_threads(Sys.CPU_THREADS)

        args = []
        trans = []
        meas = []
        for j = 1:length(N)
            push!(args, (ψtest,j))
            push!(args, (ψtest,j))
            push!(args, (ψtest,j))
            push!(args, (ψtest,j))
            push!(trans, plan_fft)
            push!(trans, plan_fft!)
            push!(trans, plan_ifft)
            push!(trans, plan_ifft!)
            push!(meas, DX[j])
            push!(meas, DX[j])
            push!(meas, DK[j])
            push!(meas, DK[j])
        end

        t1 = definetransforms(trans,args,meas,flags)

        # M = Transforms(t1...)

    return t1
end

t1 = makeTMixed(X,K)

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
