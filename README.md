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
Note
- `FourierGPE` uses `FFTW#master` built with `MKL`, pinned to a recent release. For now this _should_ avoid conflicts with `FFTW` releases, and avoid breaking packages like `FastTransforms` that require non `MKL` `FFT`.
- If you encounter an error during installation, try `]build FourierGPE`

## Getting started
To get started see runnable examples in the `/examples` directory, or for more information see [FGPEexamples.jl](https://github.com/AshtonSBradley/FGPEexamples.jl)

<img src="https://github.com/AshtonSBradley/FGPEexamples.jl/blob/master/media/3dquenchslab.gif" width="900">
