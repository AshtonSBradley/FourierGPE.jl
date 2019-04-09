# user defined paramters. Add anything here that
# you want available for GPE evaluation
@with_kw mutable struct Params @deftype Float64
    κ = 0.1
    τ = 0.3
    X::Array{Float64,2} = randn(2,2)
end

@with_kw mutable struct Sim @deftype Float64
    Lx = 200.0
    Nx::Int64 = 256
    Ly = 200.0
    Ny::Int64 = 256
    μ = 15.0
    g = 0.1
    γ = 0.5; @assert γ >= 0.0
    ti = 0.0
    tf = 1/γ
    Nt::Int64 = 200
    t::LinRange{Float64} = LinRange(ti,tf,Nt)
    ϕi::Array{Complex{Float64},2} = zeros(Nx,Ny) |> complex
    params::Params = Params() # optional parameters
    T::Transforms = Transforms()
    x::Array{Float64,1} = LinRange(0,10,100)
    y::Array{Float64,1} = LinRange(0,10,100)
    kx::Array{Float64,1} = LinRange(0,10,100)
    ky::Array{Float64,1} = LinRange(0,10,100)
    k2::Array{Complex{Float64},2} = randn(2,2) + im*randn(2,2)
end

function initsim!(sim)
# initialize transforms and fields
@unpack Lx,Ly,Nx,Ny = sim
x,y,kx,ky,k2,dx,dy,dkx,dky,Dx,Dy,Dkx,Dky,Txk,Txk!,Tkx,Tkx! = maketransforms(Lx,Ly,Nx,Ny)
T = Transforms(Txk,Txk!,Tkx,Tkx!)
@pack! sim = T,x,y,kx,ky,k2
return nothing
end

FGPETypes = [Transforms Params Sim]
for type in FGPETypes
    Base.show(io::IO, ::MIME"application/prs.juno.inline",x::type) = x
end
