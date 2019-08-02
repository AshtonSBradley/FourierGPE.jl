using Plots, SpecialFunctions

N = 500
xm = 5.
x = LinRange(-xm,xm,N); y = x

psi = @. 10*exp(-.5*(x^2+y'^2)) + randn()+im*randn()

heatmap(x,y,abs2.(psi))

# plot some binomal coefficients C_n^k
n = 1:50
plot(n,binomial.(50,n),m=:circle)

floor(10/2) |> Integer

function Gnorm(p,q,n)
    @assert q ∈ -p:2:p
    @assert n ∈ 0:2:p
    a = (p+abs(q))/2 |> Integer
    b = (p-abs(q))/2 |> Integer
    g = lfactorial(n)*lfactorial(p-n)/p/log(2)/lfactorial(a)/lfactorial(b)
    return exp(g/2)
end

Gnorm(20,4,8)

# transform function https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2858241/
function G(p,q,n)
    @assert q ∈ -p:2:p
    @assert n ∈ 0:2:p
    a = (p+abs(q))/2 |> Integer
    b = (p-abs(q))/2 |> Integer
    N = floor(n/2) |> Integer
    s = 0.0
    for r ∈ 0:N
        s += exp(lbinomial(b,r)*lbinomial(abs(q),n-2r))*im^(2r+abs(q)-n)*(-sign(q))^q
    end
    return s
end

p = 10
s = 0.
for n in 0:2:p
    global s+= G(p,-p,n)
end

s
