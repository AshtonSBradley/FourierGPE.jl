module FourierGPE

using Reexport
@reexport using OrdinaryDiffEq
@reexport using FFTW
@reexport using Parameters
@reexport using JLD2

include("arrays.jl")
include("transforms.jl")
include("evolution.jl")
include("helpers.jl")

export linspace, xvecs, kvecs, dfft
export maketransforms, makearrays, xspace, xspace!, kspace, kspace!
export nlin, nlin!, Lgp, Lgp!, V, initsim!, runsim, internalnorm
export Transforms, @pack!, @unpack, @pack_Transforms, @pack_Transforms!, @unpack_Transforms
export Par, @pack_Par, @unpack_Par

end # module
