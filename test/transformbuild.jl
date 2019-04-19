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
