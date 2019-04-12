module FourierGPE

using Reexport
@reexport using OrdinaryDiffEq
@reexport using FFTW
@reexport using Parameters
@reexport using JLD2
@reexport using FileIO

include("types.jl")
include("arrays.jl")
include("transforms.jl")
include("evolution.jl")


export Simulation, TransformLibrary, UserParams
export linspace, xvecs, kvecs, dfft, crandn_array
export maketransforms, makearrays, xspace, xspace!, kspace, kspace!
export nlin, nlin!, Lgp, Lgp!, V, initsim!, runsim, internalnorm
export Transforms, @pack_Transforms!, @unpack_Transforms
export Sim, @pack_Sim!, @unpack_Sim
export Params, @pack_Params!, @unpack_Params
export k2, @pack!, @unpack, makeT
end # module
