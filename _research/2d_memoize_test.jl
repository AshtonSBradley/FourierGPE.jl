using VortexDistributions, FourierGPE, Plots, LaTeXStrings, SpecialFunctions

using Memoize, BenchmarkTools

## CONCLUSION
# tests indicate it is still faster to call besselj0, than to memoize
# what about just preallocating?

function bessel_array(x)
    return besselj0(x)
end

## test array
L = 10.
N = 512
x = LinRange(-L,L,N)
xy = x.*x'
@btime bessel_array.(xy)

## Memoize
@memoize Dict function besselm(x)
	println("Running")
	besselj0(x)
end

# initialize
@time besselm.(xy)

# compare
@time besselm.(xy)
@time bessel_array.(xy)

## Caching
using Caching, InteractiveUtils, Serialization

@cache function besselc(x)::Float64
    besselj0(x)
end

## initialize
@time c = besselc.(xy)

## call 
@btime besselc.(xy)
