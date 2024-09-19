# FourierGPE.jl
__This project is no longer actively maintained__

[![Build Status](https://github.com/AshtonSBradley/FourierGPE.jl/workflows/CI/badge.svg)](https://github.com/AshtonSBradley/FourierGPE.jl/actions)
[![Coverage](https://codecov.io/gh/AshtonSBradley/FourierGPE.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/AshtonSBradley/FourierGPE.jl)

Simple GPE solver using FFTW. This package is intended to provide a julia implementation of the Gross-Pitaevskii equation based on Fourier spectral methods and [DifferentialEquations.jl](https://github.com/SciML/DifferentialEquations.jl)

Aims
- Simple, fast, and flexible modelling of Bose-Einstein condensate experiments.
- Allocation free, adaptive evolution in kspace.
- Establish some useful defaults and runnable examples for time evolution using `DifferentialEquations.jl`
- Arbitrary spatial dimensions.
- Not de-aliased (projective methods available elsewhere)

## Installing

```julia
]add https://github.com/AshtonSBradley/FourierGPE.jl.git
using FourierGPE
```
By default FourierGPE is built against MKL.

## Getting started
To get started see runnable examples in the `/examples` directory, or for more information see [FGPEexamples.jl](https://github.com/AshtonSBradley/FGPEexamples.jl)

<img src="https://github.com/AshtonSBradley/FGPEexamples.jl/blob/master/media/3dquenchslab.gif" width="900">
