# FourierGPE.jl
Simple GPE solver using FFTW

- Intended for fast and flexible modelling of Bose-Einstein condensate experiments
- Allocation free adaptive evolution in kspace
- Uses OrdinaryDiffEq.jl for time evolution
- Not de-aliased (projective methods available elsewhere)
- For more info see the runnable examples 

## Installing

`]add https://github.com/AshtonSBradley/FourierGPE.jl.git`

using FourierGPE

## Modifying the potential
Write your potential as a scalar function so that it can be broadcast via a dot-call, e.g.

```julia
import FourierGPE.V
V(x,y,t) = (x^2 + y^2))*pi+cos(0.1*t)
```

The potential will be called as `V.(x,y',t)`, using the column vectors `x,y` and any other parameters you add to `sim`.

## Adding parameters
You can add additional parameters (you should supply a default) to `Par` defined in 

`/src/helpers.jl`

All modified parameters will be correctly saved using, e.g. 

```julia
sim = Par()
a = 0.1
initsim!(sim)
```
which should be called prior to evolving your specified initial condition `psi` in `kspace` using

```julia
phi = kspace(psi,sim)
sol = runsim(sim,phi)
```

## Default solver
Currently uses `alg=Tsit5()`, an adaptive RK routine. The default `reltol` is not quite small enough for some applications, so it is set to `reltol = 1e-7`.
