module FourierGPE

using Reexport
@reexport using DifferentialEquations
@reexport using RecursiveArrayTools
@reexport using FFTW
@reexport using Parameters
@reexport using JLD2
@reexport using FileIO
@reexport using Plots
@reexport using LaTeXStrings

include("types.jl")
include("arrays.jl")
include("transforms.jl")
include("evolution.jl")


export Simulation, TransformLibrary, UserParams
export xvec, kvec, xvecs, kvecs, dfft, dfftall, crandn_array, crandnpartition
export maketransarrays, makearrays, xspace, xspace!, kspace, kspace!
export nlin, nlin!, Lgp, Lgp!, V, initsim!, runsim, internalnorm
export Transforms, @pack_Transforms!, @unpack_Transforms
export Sim, @pack_Sim!, @unpack_Sim, showpsi, testsim
export Params, @pack_Params!, @unpack_Params
export k2, @pack!, @unpack, makeT, makeTMixed, definetransforms
end # module
