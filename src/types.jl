
abstract type Simulation{D} end
abstract type UserParams end
abstract type TransformLibrary end

@with_kw mutable struct Transforms{D,N} <: TransformLibrary
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

@with_kw mutable struct Sim{D} <: Simulation{D} @deftype Float64
    # Add more parameters as necessary, or add to params (see examples)
    L::NTuple{D,Float64} # box length scales
    N::NTuple{D,Int64}  # grid points in each dimensions
    μ = 15.0    # chemical potential
    g = 0.1     # interaction parameter
    γ = 0.5; @assert γ >= 0.0 # damping parameter
    ti = 0.0    # initial time
    tf = 2/γ    # final time
    Nt::Int64 = 200     # number of saves over (ti,tf)
    params::UserParams # optional user parameters
    t::LinRange{Float64} = LinRange(ti,tf,Nt) # time of saves
    ϕi::Array{Complex{Float64},D} = zeros(N...) |> complex # initial condition
    alg::OrdinaryDiffEq.OrdinaryDiffEqAdaptiveAlgorithm = Tsit5() # default solver
    reltol::Float64 = 1e-6 # default tolerance; may need to use 1e-7 for corner cases
    flags::UInt32 = FFTW.MEASURE # choose a plan. PATIENT, NO_TIMELIMIT, EXHAUSTIVE
    nfiles::Bool = false
    path::String = nfiles ? "data" : "./"
    filename::String = "save"
    # =======================================
    # arrays, transforms, spectral operators
    X::NTuple{D,Array{Float64,1}} = xvecs(L,N)
    K::NTuple{D,Array{Float64,1}} = kvecs(L,N)
    espec::Array{Complex{Float64},D} = 0.5*k2(K)
    T::TransformLibrary = makeT(X,K,flags=flags)
end

Sim(L,N,par) = Sim{length(L)}(L=L,N=N,params=par)

FGPETypes = [Transforms Sim]
for type in FGPETypes
    Base.show(io::IO, ::MIME"application/prs.juno.inline",x::type) = x
end
