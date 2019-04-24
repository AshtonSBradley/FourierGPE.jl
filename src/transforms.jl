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
    T = makeT(X,K,N)

Build transform library for the array tuples `X`, `K`. Defaults to a measure plan.
`j` is number of scratch fields to initialize for in-place evaluation.
"""
function makeT(X,K,j=1;flags=FFTW.MEASURE)
    FFTW.set_num_threads(Sys.CPU_THREADS)
    dim = length(X)
    N = length.(X)
    DX,DK = dfftall(X,K)
    dμx = prod(DX)
    dμk = prod(DK)
    ψtest = ones(N...) |> complex

    trans = (plan_fft,plan_fft!,plan_ifft,plan_ifft!)
    meas = (dμx,dμx,dμk,dμk)
    flags = FFTW.MEASURE
    args = ((ψtest,),(ψtest,),(ψtest,),(ψtest,))
    Txk,Txk!,Tkx,Tkx! = definetransforms(trans,args,meas,flags)
    Mxk,Mxk!,Mkx,Mkx! = makeTMixed(X,K,flags=flags)

    return Transforms{dim,j}(Txk,Txk!,Tkx,Tkx!,Mxk,Mxk!,Mkx,Mkx!,crandnpartition(dim,j))
end

function makeTMixed(X,K;flags=FFTW.MEASURE)
    FFTW.set_num_threads(Sys.CPU_THREADS)
    N = length.(X)
    dim = length(X)
    DX,DK = dfftall(X,K)
    ψtest = ones(N...) |> complex

args = []
transxk = []
transxk! = []
transkx = []
transkx! = []
measxk = []
measkx = []


for i = 1:dim
    push!(args, (ψtest,i))
    push!(transxk, plan_fft)
    push!(transxk!, plan_fft!)
    push!(transkx, plan_ifft)
    push!(transkx!, plan_ifft!)
    push!(measxk, DX[i])
    push!(measkx, DK[i])

end
Mxk = definetransforms(transxk,args,measxk,flags)
Mxk! = definetransforms(transxk!,args,measxk,flags)
Mkx = definetransforms(transkx,args,measkx,flags)
Mkx! = definetransforms(transkx!,args,measkx,flags)
return Mxk,Mxk!,Mkx,Mkx!
end

# function makeT(X,K,flags=FFTW.MEASURE)
#
#     N = length.(X)
#     DX,DK = dfftall(X,K)
#     dμx = prod(DX)
#     dμk = prod(DK)
#     ψtest = ones(N...) |> complex
#
#     trans = (plan_fft,plan_fft!,plan_ifft,plan_ifft!)
#     meas = (dμx,dμx,dμk,dμk)
#     flags = FFTW.MEASURE
#     args = ((ψtest,),(ψtest,),(ψtest,),(ψtest,))
#     FFTW.set_num_threads(Sys.CPU_THREADS)
#
#     return definetransforms(trans,args,meas,flags)
# end

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

function maketransarrays(L,N,j=1;flags=FFTW.MEASURE)

    X,K,dX,dK = makearrays(L,N)
    T = makeT(X,K,j,flags=flags)

return X,K,dX,dK,DX,DK,T
end
