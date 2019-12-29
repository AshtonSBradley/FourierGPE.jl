"""
    Dx,Dk = dfft(x,k)

Measures that make `fft`, `ifft` 2-norm preserving.
Correct measures for mapping between `x`- and `k`-space.
"""
function dfft(x,k)
    dx = diff(x)[1]; dk = diff(k)[1]
    Dx = dx/sqrt(2*pi)
    Dk = length(k)*dk/sqrt(2*pi)
    return Dx, Dk
end

"""
    DX,DK = dfftall(X,K)

Evalutes tuple of measures that make `fft`, `ifft` 2-norm preserving for each
`x` or `k` dimension.
"""
function dfftall(X,K)
    M = length(X)
    DX = zeros(M); DK = zeros(M)
    for i ∈ eachindex(X)
        DX[i],DK[i] = dfft(X[i],K[i])
    end
    return DX,DK
end

"""
    ψ = xspace(ϕ,sim)

Transform from `k`- to `x`-space using transforms packed into `sim`.
"""
function xspace(ϕ,sim)
    @unpack T = sim
    return T.Tkx*ϕ
end

"""
    xspace!(ϕ,sim)

Mutating transform from `k`- to `x`-space using transforms packed into `sim`.
"""
function xspace!(ψ,sim)
    @unpack T = sim
    T.Tkx!*ψ
    return nothing
end

"""
    kspace(ψ,sim)

Transform from `x`- to `k`-space using transforms packed into `sim`.
"""
function kspace(ψ,sim)
    @unpack T = sim
    return T.Txk*ψ
end

"""
    kspace!(ψ,sim)

Mutating transform from `x`- to `k`-space using transforms packed into `sim`.
"""
function kspace!(ψ,sim)
    @unpack T = sim
    T.Txk!*ψ
    return nothing
end

"""
    definetransforms(funcs,args,meas,kwargs)

Build all transforms for simulation.
"""
function definetransforms(funcs,args,meas,kwargs)
    trans = []
    for (fun,arg) in zip(funcs,args)
        push!(trans, fun(arg...,flags=kwargs))
    end
    return meas.*trans
end

"""
    T = makeT(X,K,j)

Build `FFTW` transform library for the array tuples `X`, `K`. Defaults to a measure plan.
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
    # flags = FFTW.MEASURE
    args = ((ψtest,),(ψtest,),(ψtest,),(ψtest,))

    # flags != FFTW.MEASURE && @info "Planning FFTs ..."
    Txk,Txk!,Tkx,Tkx! = definetransforms(trans,args,meas,flags)
    # Mxk,Mxk!,Mkx,Mkx! = makeTMixed(X,K,flags=flags)
    # @info "...Plans created."
    return Transforms{dim,j}(Txk,Txk!,Tkx,Tkx!,crandnpartition(dim,j))
end

"""
    T = makeTMixed(X,K,j)

Build mixed transform library for the array tuples `X`, `K`. Defaults to a measure plan.
"""
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

#TODO fails in MKL
Mxk = definetransforms(transxk,args,measxk,flags)
Mxk! = definetransforms(transxk!,args,measxk,flags)
Mkx = definetransforms(transkx,args,measkx,flags)
Mkx! = definetransforms(transkx!,args,measkx,flags)
return Mxk,Mxk!,Mkx,Mkx!
end

"""
    maketransarrays(L,N,j=1;flags=FFTW.MEASURE)

Make all transforms and arrays in one call.
"""
function maketransarrays(L,N,j=1;flags=FFTW.MEASURE)
    X,K,dX,dK = makearrays(L,N)
    T = makeT(X,K,j,flags=flags)
    return X,K,dX,dK,DX,DK,T
end
