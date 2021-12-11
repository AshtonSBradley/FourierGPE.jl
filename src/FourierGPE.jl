module FourierGPE

using Reexport
using SpecialFunctions
using LazyArrays
using FillArrays
using OrdinaryDiffEq
using LoopVectorization    
using Tullio               
using RecursiveArrayTools
using PaddedViews
using FFTW
import DiffEqCallbacks.FunctionCallingCallback

# using LaTeXStrings
# using ColorSchemes
# using Plots
# const c1 = cgrad(ColorSchemes.inferno.colors)
# const c2 = cgrad(ColorSchemes.RdBu_11.colors)
# const c3 =:mediumseagreen

@reexport using Parameters
@reexport using JLD2
@reexport using FileIO

include("types.jl")
include("arrays.jl")
include("transforms.jl")
include("evolution.jl")

# simulation
export Simulation, TransformLibrary, UserParams
export xvec, kvec, xvecs, kvecs, dfft, dfftall, crandn_array, crandnpartition 
export maketransarrays, makearrays, xspace, xspace!, kspace, kspace!
export nlin, nlin!, Lgp, Lgp!, V, Params
export initsim!, runsim, internalnorm
export Transforms, @pack_Transforms!, @unpack_Transforms
export Sim, @pack_Sim!, @unpack_Sim, testsim
# export showpsi, c1, c2, c3
export Params, @pack_Params!, @unpack_Params
export k2, @pack!, @unpack, makeT, definetransforms #makeTMixed,
export Field, XField, KField

end # module
