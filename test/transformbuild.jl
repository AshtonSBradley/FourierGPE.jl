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
    for j ∈ eachindex(funcs)
        push!(trans, funcs[j](args...,flags=kwargs))
    end
    return meas.*trans
end

t1 = deftrans(trans,args,flags)


# ====== simpler approach =====

T = Transforms(t1...)


# Txk = dμx*plan_fft(ψtest,flags=flags)
# # ψtest = ones(Ns...) |> complex
# Txk! = dμx*plan_fft!(ψtest,flags=flags)
# # ψtest = ones(Ns...) |> complex
# Tkx  = dμk*plan_ifft(ψtest,flags=flags)
# # ψtest = ones(Ns...) |> complex
# Tkx!  = dμk*plan_ifft!(ψtest,flags=flags)
# # ψtest = ones(Ns...) |> complex
