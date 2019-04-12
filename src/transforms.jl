# transform methods

"""Measures that make fft, ifft 2-norm preserving."""
function dfft(x,k)
    dx = diff(x)[1]; dk = diff(k)[1]
    Dx = dx/sqrt(2*pi)
    Dk = length(k)*dk/sqrt(2*pi)
    return Dx, Dk
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

function maketransforms(L,N;flags=FFTW.MEASURE)
    X = xvecs(L,N)
    K = kvecs(L,N)
    M = length(N)
    #measures for unitary FFT, and for standard integration
    dX = zeros(M); dK = zeros(M)
    DX = zeros(M); DK = zeros(M)
    for i in eachindex(X)
        DX[i],DK[i] = dfft(X[i],K[i])
        dX[i],dK[i] = diff(X[i])[1],diff(K[i])[1]
    end
    dμx = prod(DX); dμk = prod(DK)

    # plan transforms
    FFTW.set_num_threads(Sys.CPU_THREADS)
    ψtest = ones(N...) |> complex
    Txk = dμx*plan_fft(ψtest,flags=flags)
    ψtest = ones(N...) |> complex
    Txk! = dμx*plan_fft!(ψtest,flags=flags)
    ψtest = ones(N...) |> complex
    Tkx  = dμk*plan_ifft(ψtest,flags=flags)
    ψtest = ones(N...) |> complex
    Tkx!  = dμk*plan_ifft!(ψtest,flags=flags)
    ψtest = ones(N...) |> complex
    T = Transforms(Txk,Txk!,Tkx,Tkx!)
return X,K,dX,dK,DX,DK,T
end

"""
    T = makeT(X,K;flags=FFTW.MEASURE)

Build transform library for the array tuples `X`, `K`. Defaults to a measure plan.
"""
function makeT(X,K;flags=FFTW.MEASURE)
    M = length(X)
    N = [ length(X[i]) for i ∈ eachindex(X) ] |> Tuple
    #measures for unitary FFT
    dX = zeros(M); dK = zeros(M)
    DX = zeros(M); DK = zeros(M)
    for i ∈ eachindex(X)
        DX[i],DK[i] = dfft(X[i],K[i])
        dX[i],dK[i] = diff(X[i])[1],diff(K[i])[1]
    end
    dμx = prod(DX); dμk = prod(DK)

    # plan transforms
    FFTW.set_num_threads(Sys.CPU_THREADS)
    ψtest = ones(N...) |> complex
    Txk = dμx*plan_fft(ψtest,flags=flags)
    ψtest = ones(N...) |> complex
    Txk! = dμx*plan_fft!(ψtest,flags=flags)
    ψtest = ones(N...) |> complex
    Tkx  = dμk*plan_ifft(ψtest,flags=flags)
    ψtest = ones(N...) |> complex
    Tkx!  = dμk*plan_ifft!(ψtest,flags=flags)
    ψtest = ones(N...) |> complex
    T = Transforms(Txk,Txk!,Tkx,Tkx!)
    return T
end
