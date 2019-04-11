
abstract type Simulation{D} end
abstract type UserParams end
abstract type TransformLibrary{D} end

@with_kw mutable struct Transforms{D} <: TransformLibrary{D}
    Txk::AbstractFFTs.ScaledPlan{Complex{Float64},FFTW.cFFTWPlan{Complex{Float64},-1,false,D},Float64} = 0.1*plan_fft(crandn_array(D))
    Txk!::AbstractFFTs.ScaledPlan{Complex{Float64},FFTW.cFFTWPlan{Complex{Float64},-1,true,D},Float64} = 0.1*plan_fft!(crandn_array(D))
    Tkx::AbstractFFTs.ScaledPlan{Complex{Float64},FFTW.cFFTWPlan{Complex{Float64},1,false,D},Float64} = 0.1*plan_ifft(crandn_array(D))
    Tkx!::AbstractFFTs.ScaledPlan{Complex{Float64},FFTW.cFFTWPlan{Complex{Float64},1,true,D},Float64} = 0.1*plan_ifft!(crandn_array(D))
end

@with_kw mutable struct Sim{D} <: Simulation{D} @deftype Float64
    L::NTuple{D,Float64}
    N::NTuple{D,Int64}
    μ = 15.0
    g = 0.1
    γ = 0.5; @assert γ >= 0.0
    ti = 0.0
    tf = 2/γ
    Nt::Int64 = 200
    t::LinRange{Float64} = LinRange(ti,tf,Nt)
    ϕi::Array{Complex{Float64},D} = zeros(N...) |> complex
    params::UserParams # optional parameters
    X::NTuple{D,Array{Float64,1}} = xvecs(L...,N...)
    K::NTuple{D,Array{Float64,1}} = kvecs(L...,N...)
    espec::Array{Complex{Float64},D} = k2(L...,N...)
    T::TransformLibrary{D} = Transforms{D}()
end

Sim(L,N,par) = Sim{length(L)}(L=L,N=N,params=par)

FGPETypes = [Transforms Sim]
for type in FGPETypes
    Base.show(io::IO, ::MIME"application/prs.juno.inline",x::type) = x
end
