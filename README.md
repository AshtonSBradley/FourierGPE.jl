# FourierGPE.jl
Simple GPE solver using FFTW

- Intended to provide simple, fast, and flexible modelling of Bose-Einstein condensate experiments.
- Allocation free, adaptive evolution in kspace.
- Establish some useful defaults and runnable examples for time evolution using `DifferentialEquations.jl`
- Arbitrary spatial dimensions.
- Not de-aliased (projective methods available elsewhere)

## Installing

```julia
]add https://github.com/AshtonSBradley/FourierGPE.jl.git
using FourierGPE
```
There is an option to build against MKL, but this will cause issues with 3D simulations. Hence, if you want MKL enhancements for 1D,2D, remove the comment  in `build.jl`, and rebuild. 

## Getting started
To get started see runnable examples in the `/examples` directory, or for more information see [FGPEexamples.jl](https://github.com/AshtonSBradley/FGPEexamples.jl)

<img src="https://github.com/AshtonSBradley/FGPEexamples.jl/blob/master/media/3dquenchslab.gif" width="900">
