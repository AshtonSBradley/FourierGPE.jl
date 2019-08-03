using Test, Plots, Parameters

Nx = 500
xm = 10.
x = LinRange(-xm,xm,Nx); y = x

# create some test data
psi = @. 10*exp(-.1*(x^2+y'^2)) + 15*exp(-((x-5)^2+y'^2)) # + randn()+im*randn()
heatmap(x,y,abs2.(psi))

abstract type Basis end
struct Oscillator <: Basis
    ω::Float64
    n::Int64
end
Oscillator(n::Int64) = Oscillator(1.0,n)
function (b::Oscillator)(x)
    @unpack ω,n = b
    T = x * ones(1, n) |> zero
    T[:,1] = @. exp(-(√ω * x)^2 / 2) * (ω / π)^(1 / 4)
    n >= 2 && (T[:,2] = @. sqrt(2 * ω) * x * T[:,1])
    for j ∈ index(b)[3:end]
        m = qnumbers(b)[j]
        T[:,j] = @. sqrt(2 / m) * (√ω * x) * T[:,j-1] - sqrt((m-1) / m) * T[:,j-2]
    end
    return T
end
index(b::Oscillator) = 1:b.n
spectrum(b::Oscillator) = (0:b.n-1)*b.ω
qnumbers(b::Oscillator) = 0:b.n-1

function hermite(x,n)
    @assert n >= 1
    h1 = exp(-x^2 / 2) * (1 / π)^(1 / 4)
    n == 1 && (return h1)
    h2 = sqrt(2) * x * h1
    n == 2 && (return h2)
    h3 = 0.0
    for j ∈ 3:n
        m = j - 1
        h3 = sqrt(2/m) * x * h2 - sqrt((m-1)/m) * h1
        h1 = h2; h2 = h3
    end
    n > 2 && (return h3)
end
hermite(x,n,ω) = hermite(√ω * x,n) * ω ^ ( 1 / 4 )

# test that we can create the basis, or any one of oscillator modes
Nv = 1:20
@time hermite.(x,Nv')
plot(x,hermite.(x,Nv'),legend=false,grid=false)
n = 20
b1 = Oscillator(n)
@time b1(x)
plot(x,b1(x),legend=false,grid=false)

plot(x,hermite.(x,1,[1. .5 .2 .1]))

# filtering cartesian data
function filterH(psi,x,N)
    H = Oscillator(1.0,N)(x)
    T = inv(H'*H)
    F̂ = T*H'*psi*H*T
    psiF = zero(psi)
        for m = 1:N, n = 1:(N+1-m)
            @. psiF += H[:,m]*F̂[m,n]*H[:,n]'
        end
    return psiF
end

n = 30
psiF = filterH(psi,x,n)
heatmap(x,y,abs2.(psiF))

# coversion to polar coords
function polar(psi,x,N)
    Nx = length(x)
    H = Oscillator(1.0,N)(x)
    T = inv(H'*H)
    F̂ = T*H'*psi*H*T
    r = LinRange(0,last(x),Nx/2 |> Int)
    θ = LinRange(0,2*pi,2*Nx)'
    psiP = zero(r*θ)
    for m = 1:N, n = 1:(N + 1 - m)
        @. psiP += hermite(r*cos(θ),n) * F̂[m,n] * hermite(r*sin(θ),m)
    end
    return psiP
end

r = LinRange(0,last(x),Nx/2 |> Int)
θ = LinRange(0,2*pi,2*Nx)'
psiP = polar(psi,x,30)

# plot in polar coordinates
heatmap(θ',r,abs2.(psiP))
dθ = θ[2]-θ[1]
dr = r[2]-r[1]
dx = x[2]-x[1]
psiPth = sum(psiP,dims=2)*dθ
plot(r,psiPth.*r)

# check total integral preserved (< 1%)
sum(psiPth.*r)*dr
sum(psiF)*dx^2

# ======== ideas for cfields
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
