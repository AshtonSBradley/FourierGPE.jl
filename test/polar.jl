using Test, Plots, SpecialFunctions, FastGaussQuadrature, Parameters
gr(aspect_ratio=1)

Nx = 500
xm = 10.
x = LinRange(-xm,xm,Nx); y = x
psi = @. 10*exp(-.5*(x^2+y'^2)) + randn()+im*randn()
heatmap(x,y,abs2.(psi))

abstract type Basis end
struct Oscillator <: Basis
    ω::Float64
    n::Int64
end
Oscillator(n::Int64) = Oscillator(1.0,n)

b1 = Oscillator(10)

index(b::Oscillator) = 1:b.n
spectrum(b::Oscillator) = (0:b.n-1)*b.ω
qnumbers(b::Oscillator) = 0:b.n-1

function (b::Oscillator)(x)
    @unpack ω,n = b
    T = x * ones(1, n) |> zero
    T[:,1] = @. exp(-(√ω * x)^2 / 2) * (ω / π)^(1 / 4)
    n >= 2 && (T[:,2] = @. sqrt(2) * exp(-(√ω * x)^2 / 2) * (√ω * x) * (ω / π)^(1 / 4))
    for j ∈ index(b)[3:end]
        m = qnumbers(b)[j]
        T[:,j] = @. sqrt(2 / m) * (√ω * x) * T[:,j-1] - sqrt((m-1) / m) * T[:,j-2]
    end
    return T
end



n = 12
b1 = Oscillator(1.,n)
plot(x,b1(x),legend=false,grid=false)



# CField:

#= something like this?
function transform(b::Oscillator)
    x,w = gausshermite(b.n)
    T = b(x)
    return x,w,T
end

xk,wk,Tk = transform(b1)
c = randn(n)
f = Tk*c
@test sum(@. wk*exp(xk^2)*abs2(f)) ≈ sum(abs2.(c))

abstract type CField end
struct YField{N} <: CField
    basis::Basis
    psi::Array{Complex{Float64},N}
end

f1 = YField(b1,randn(b1.n)+im*randn(b1.n))

=#
