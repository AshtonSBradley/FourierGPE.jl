module FourierGPE

using Reexport
using SpecialFunctions
using LazyArrays
using FillArrays
using DifferentialEquations
using RecursiveArrayTools
using FFTW
using LaTeXStrings
using ColorSchemes
@reexport using Parameters
@reexport using JLD2
@reexport using FileIO
import Plots
const c1 = Plots.cgrad(ColorSchemes.inferno.colors)
const c2 = Plots.cgrad(ColorSchemes.RdBu_11.colors)
const c3 =:mediumseagreen

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
export Sim, @pack_Sim!, @unpack_Sim, testsim
export showpsi, c1, c2, c3
export Params, @pack_Params!, @unpack_Params
export k2, @pack!, @unpack, makeT, definetransforms #makeTMixed,
export Field, XField, KField, velocity, energydecomp, helmholtz
export log10range, zeropad, autocorrelate, convolve, kespectrum, ikespectrum
end # module
