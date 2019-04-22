# transform methods

"""Measures that make fft, ifft 2-norm preserving."""
function dfft(x,k)
    dx = diff(x)[1]; dk = diff(k)[1]
    Dx = dx/sqrt(2*pi)
    Dk = length(k)*dk/sqrt(2*pi)
    return Dx, Dk
end

function dfftall(X,K)
    M = length(X)
    DX = zeros(M); DK = zeros(M)
    for i ∈ eachindex(X)
        DX[i],DK[i] = dfft(X[i],K[i])
    end
    return DX,DK
end

    function xspace(ψ,sim)
        @unpack T = sim
        return T.Tkx*ψ
    end

    function xspace!(ψ,sim)
        @unpack T = sim
        T.Tkx!*ψ
        return nothing
    end

    function kspace(ψ,sim)
        @unpack T = sim
        return T.Txk*ψ
    end

    function kspace!(ψ,sim)
        @unpack T = sim
        T.Txk!*ψ
        return nothing
    end

function definetransforms(funcs,args,meas,kwargs)
    trans = []
    for (i,fun) ∈ enumerate(funcs)
        push!(trans, fun(args[i]...,flags=kwargs))
    end
    return meas.*trans
end

"""
    T = makeT(X,K)

Build transform library for the array tuples `X`, `K`. Defaults to a measure plan.
"""
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

# function makeT(X,K,flags=FFTW.MEASURE)
#
#     N = [ length(X[i]) for i ∈ eachindex(X) ] |> Tuple
#     DX,DK = dfftall(X,K)
#     dμx = prod(DX)
#     dμk = prod(DK)
#     ψtest = ones(N...) |> complex
#
#     # plan transforms
#     FFTW.set_num_threads(Sys.CPU_THREADS)
#     Txk = dμx*plan_fft(ψtest,flags=flags)
#     Txk! = dμx*plan_fft!(ψtest,flags=flags)
#     Tkx  = dμk*plan_ifft(ψtest,flags=flags)
#     Tkx!  = dμk*plan_ifft!(ψtest,flags=flags)
#     T = Transforms(Txk,Txk!,Tkx,Tkx!)
#
#     return T
# end

function maketransarrays(L,N,flags=FFTW.MEASURE)

    X,K,dX,dK = makearrays(L,N)
    T = Transforms(makeT(X,K,flags)...)

return X,K,dX,dK,DX,DK,T
end
