module FourierGPE

using Reexport
@reexport using SpecialFunctions
@reexport using LazyArrays
@reexport using FillArrays
@reexport using DifferentialEquations
@reexport using RecursiveArrayTools
@reexport using FFTW
@reexport using Parameters
@reexport using JLD2
@reexport using FileIO
@reexport using Plots
@reexport using LaTeXStrings
@reexport using ColorSchemes

const c2 = cgrad(ColorSchemes.RdBu_11.colors)
const c1 = cgrad(ColorSchemes.linear_blue_5_95_c73_n256.colors)

include("types.jl")
include("arrays.jl")
include("transforms.jl")
include("evolution.jl")
include("analysis.jl")

export Simulation, TransformLibrary, UserParams
export xvec, kvec, xvecs, kvecs, dfft, dfftall, crandn_array, crandnpartition
export maketransarrays, makearrays, xspace, xspace!, kspace, kspace!
export nlin, nlin!, Lgp, Lgp!, V, Params
export initsim!, runsim, internalnorm
export Transforms, @pack_Transforms!, @unpack_Transforms
export Sim, @pack_Sim!, @unpack_Sim, showpsi, testsim
export Params, @pack_Params!, @unpack_Params
export k2, @pack!, @unpack, makeT, makeTMixed, definetransforms
export Field, XField, KField, velocity, energydecomp, helmholtz
export log10range, zeropad, autocorrelate, convolve, kespectrum, ikespectrum
end # module
