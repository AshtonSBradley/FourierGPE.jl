
abstract type Simulation{D} end
abstract type UserParams end
abstract type TransformLibrary end
abstract type Field end

struct XField{D} <: Field
    psiX::Array{Complex{Float64},D}
    X::NTuple{D}
    K::NTuple{D}
    K2::Array{Float64,D}
end

struct KField{D} <: Field
    psiK::Array{Complex{Float64},D}
    X::NTuple{D}
    K::NTuple{D}
    K2::Array{Float64,D}
end

@with_kw mutable struct Transforms{D,N} <: TransformLibrary
    Txk::AbstractFFTs.ScaledPlan{Complex{Float64},FFTW.cFFTWPlan{Complex{Float64},-1,false,D,UnitRange{Int64}},Float64} = 0.1*plan_fft(crandn_array(D))
    Txk!::AbstractFFTs.ScaledPlan{Complex{Float64},FFTW.cFFTWPlan{Complex{Float64},-1,true,D,UnitRange{Int64}},Float64} = 0.1*plan_fft!(crandn_array(D))
    Tkx::AbstractFFTs.ScaledPlan{Complex{Float64},FFTW.cFFTWPlan{Complex{Float64},1,false,D,UnitRange{Int64}},Float64} = 0.1*plan_ifft(crandn_array(D))
    Tkx!::AbstractFFTs.ScaledPlan{Complex{Float64},FFTW.cFFTWPlan{Complex{Float64},1,true,D,UnitRange{Int64}},Float64} = 0.1*plan_ifft!(crandn_array(D))
    psi::ArrayPartition = crandnpartition(D,N)
end

# ==== define user parameters =======
# rename and make this whatevery you require and pass it to sim
@with_kw mutable struct Params <: UserParams @deftype Float64
    κ = 0.0 # a placeholder
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
    params::UserParams = Params() # optional user parameters
    V0::Array{Float64,D} = zeros(N)
    t::LinRange{Float64} = LinRange(ti,tf,Nt) # time of saves
    ϕi::Array{Complex{Float64},D} = zeros(N) |> complex # initial condition
    alg::OrdinaryDiffEq.OrdinaryDiffEqAdaptiveAlgorithm = Tsit5() # default solver
    reltol::Float64 = 1e-6 # default tolerance; may need to use 1e-7 for corner cases
    flags::UInt32 = FFTW.MEASURE # choose a plan. PATIENT, NO_TIMELIMIT, EXHAUSTIVE
    # === saving
    nfiles::Bool = false
    path::String = nfiles ? joinpath(@__DIR__,"data") : @__DIR__
    filename::String = "save"
    # === arrays, transforms, spectral operators
    X::NTuple{D,Array{Float64,1}} = xvecs(L,N)
    K::NTuple{D,Array{Float64,1}} = kvecs(L,N)
    espec::Array{Complex{Float64},D} = 0.5*k2(K)
    T::TransformLibrary = makeT(X,K,flags=flags)
end

Sim(L,N,par) = Sim{length(L)}(L=L,N=N,params=par)
Sim(L,N) = Sim{length(L)}(L=L,N=N,params=Params())

FGPETypes = [Transforms Sim]
for type in FGPETypes
    Base.show(io::IO, ::MIME"application/prs.juno.inline",x::type) = x
end
