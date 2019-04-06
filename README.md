# FourierGPE.jl
Simple GPE solver using FFTW

- Intended to provide simple, fast, and flexible modelling of Bose-Einstein condensate experiments
- Allocation free adaptive evolution in kspace
- Establish some useful defaults for time evolution in `OrdinaryDiffEq.jl`
- Not de-aliased (projective methods available elsewhere)
- For more info see the runnable examples 

## Installing

```julia
]add https://github.com/AshtonSBradley/FourierGPE.jl.git
using FourierGPE
```

## Modifying the potential
Write your potential as a scalar function so that it can be broadcast via a dot-call, e.g.

```julia
import FourierGPE.V
V(x,y,t) = (x^2 + y^2))*pi+cos(0.1*t)
```

The potential will be called as `V.(x,y',t)`, using the column vectors `x,y` and any other parameters you add to `sim`.

## Adding parameters
You can add additional parameters by editing the definition of `Par` in `/src/helpers.jl`,
which will require you to supply a default, and optionally a type (if it differs from `Float64`).

Modified parameters for a particular simulation will be correctly saved and made available using, e.g. 

```julia
sim = Par() # initialize default simulation
a = 0.1 # modified value 
@pack! sim = a #pack it into sim
initsim!(sim)  # create initial arrays and transforms
@unpack_Par sim # provde all variables in the current workspace
```
which should be called prior to evolving your specified initial condition `psi` in `kspace` using

```julia
phi = kspace(psi,sim)
sol = runsim(sim,phi)
```

## Default solver
Currently uses `alg=Tsit5()`, an adaptive RK routine. The default `reltol` is not quite small enough for some applications, so it is set to `reltol = 1e-7`. See `src/evolution.jl` for details.
