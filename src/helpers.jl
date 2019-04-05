import Base.randn
randn(a::Array{T,N}) where {T,N} = randn(T,size(a))

@with_kw mutable struct Par @deftype Float64
    Lx = 200.0
    Nx::Int64 = 256
    Ly = 200.0
    Ny::Int64 = 256
    μ = 15.0
    g = 0.1
    γ = 0.5
    ti = 0.0
    tf = 1/γ
    Nt::Int64 = 200
    t::LinRange{Float64} = LinRange(ti,tf,Nt)
    T::Trans = Trans()
    x::Array{Float64,1} = LinRange(0,10,100)
    y::Array{Float64,1} = LinRange(0,10,100)
    kx::Array{Float64,1} = LinRange(0,10,100)
    ky::Array{Float64,1} = LinRange(0,10,100) 
    k2::Array{Complex{Float64},2} = randn(2,2) + im*randn(2,2)
end

function initsim!(sim)
# initialize transforms and fields
@unpack Lx,Nx,Ly,Ny = sim
x,y,kx,ky,k2,dx,dy,dkx,dky,Dx,Dy,Dkx,Dky,Txk,Txk!,Tkx,Tkx! = maketransforms(Lx,Nx,Ly,Ny)
T = Trans(Txk,Txk!,Tkx,Tkx!)
@pack! sim = T,x,y,kx,ky,k2
return nothing
end
