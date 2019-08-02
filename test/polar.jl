using Test, Plots, SpecialFunctions, FastGaussQuadrature, Parameters

N = 500
xm = 5.
x = LinRange(-xm,xm,N); y = x
psi = @. 10*exp(-.5*(x^2+y'^2)) + randn()+im*randn()
heatmap(x,y,abs2.(psi))

M = 20
x,w = gausshermite(M)
Ω = 1.0

abstract type Basis end
struct Oscillator <: Basis
    Ω::Float64
    n::Int64
end
Oscillator() = Oscillator(1.0,10)
b1 = Oscillator()


# Oscillator(n) = Oscillator(1.,n)
# Oscillator() = Oscillator(10)
# osc = Oscillator()
# Oscillator(n::Vector) = Oscillator.(n)
# Oscillator(n::UnitRange{Int64}) = Oscillator.(n)

# n = 0:10

# Oscillator.((10,20))

# function eigenstatematrix(basis::Oscillator{1},x)
    # @unpack_Oscillator basis
    # M = N[1];  Ω = ω[1]
# @assert M >= 1
# M > 370 && error("Float64 Quadrature does not converge for M > 370.")


function (t::Oscillator)(x)
    @unpack Ω,n = t
    m = 0:n - 1 |> Vector
    T = x * ones(1, n) |> zero
    T[:,1] = @. exp(-(√Ω * x)^2 / 2) * (Ω / π)^(1 / 4)
    n >= 2 && (T[:,2] = @. sqrt(2) * exp(-(√Ω * x)^2 / 2) * (√Ω * x) * (Ω / π)^(1 / 4))
    for j = 1:n - 2
        T[:,j + 2] = @. sqrt(2 / (m[j + 2])) * (√Ω * x) * T[:,j + 1] - sqrt(m[j + 1] / m[j + 2]) * T[:,j]
    end
    return T
end

b1 = Oscillator(1.,20)
T = b1(x)
c = randn(20)
f = T*c

@test sum(@. w*exp(x^2)*abs2(f)) ≈ sum(abs2.(c))

#something like this?
function transform(b::Oscillator)
    x,w = gausshermite(b.n)
    T = b(x)
    return x,w,T
end
