# array methods
linspace(a,b,n) = LinRange(a,b,n) |> collect

xvec(L,N)=LinRange(-L/2,L/2,N+1)[2:end] |> collect

function kvec(λ,N)
    @assert iseven(N)
    nkx = 0:Int(N/2)
    kx = [nkx[1:end-1];-reverse(nkx[2:end])]*2*π/λ
return kx
end

function xvecs(L,N)
    X = []
    for (λ,ν) in zip(L,N)
        x = xvec(λ,ν)
        push!(X,x)
    end
    return X |> Tuple
end

function kvecs(L,N)
    K = []
    for (λ,ν) in zip(L,N)
        k = kvec(λ,ν)
        push!(K,k)
    end
    return K |> Tuple
end

function k2(K)
    kind = Iterators.product(K...)
    return map(x-> sum(abs2.(x)),kind)
end

function makearrays(L,N)
    X = xvecs(L,N)
    K = kvecs(L,N)
    dX = Float64[]; dK = Float64[]
    for j ∈ eachindex(X)
        x=X[j];k=K[j]
        dx = x[2]-x[1];dk = k[2]-k[1]
        push!(dX,dx)
        push!(dK,dk)
    end
    dX = dX |> Tuple
    dK = dK |> Tuple
    return X,K,dX,dK
end

"""
A = crandn_array(M)

Make placeholder `2x2x...` complex `randn()` array of `M` dimensions."""
function crandn_array(M)
    a = Int.(ones(M)).+1
    return randn(a...) |> complex
end

"""
A = crandnpartition(D,M)

Make placeholder ArrayPartition vector of length `M`, containing `2x2x...` rank D complex matrices.
"""
function crandnpartition(D,N)
    a = crandn_array(D)
    args = []
    for j = 1:N
        push!(args,a)
    end
    return ArrayPartition(args...)
end
