# array methods
linspace(a,b,n) = LinRange(a,b,n) |> collect

"""
    x = xvec(λ,N)

Creates `x` values with correct periodicity for box specified by length `λ` for number of points `N`.
"""
xvec(L,N) = LinRange(-L/2,L/2,N+1)[2:end] |> collect

"""
    k = kvec(λ,N)

Creates `k` values with correct periodicity for box specified by length `λ` for number of points `N`.
"""
function kvec(λ,N)
    @assert iseven(N)
    nk = 0:Int(N/2)
    k = [nk[1:end-1];-reverse(nk[2:end])]*2*π/λ
    return k
end

"""
    X = xvecs(L,N)

Creates a tuple containing the spatial coordinate array for each spatial dimension.
"""
function xvecs(L,N)
    X = []
    for (λ,ν) in zip(L,N)
        x = xvec(λ,ν)
        push!(X,x)
    end
    return X |> Tuple
end

"""
    K = kvecs(L,N)

Creates a tuple containing the spatial coordinate array for each spatial dimension.
"""
function kvecs(L,N)
    K = []
    for (λ,ν) in zip(L,N)
        k = kvec(λ,ν)
        push!(K,k)
    end
    return K |> Tuple
end

"""
    k² = k2(K)

Creates the kinetic energy array `k²` on the `k`-grid defined by the tuple `K`.
"""
function k2(K)
    kind = Iterators.product(K...)
    return map(x-> sum(abs2.(x)),kind)
end

"""
    X,K,dX,dK = makearrays(L,N)

Creates all `x` and `k` arrays for box specified by tuples `L=(Lx,...)` and `N=(Nx,...)`.
Differenetials `dX`, `dK` are also reaturned. `L` and `N` must be of equal length.
"""
function makearrays(L,N)
    @assert length(L) == length(N)
    X = xvecs(L,N)
    K = kvecs(L,N)
    dX = Float64[]; dK = Float64[]
    for j ∈ eachindex(X)
        x = X[j]; k = K[j]
        dx = x[2]-x[1]; dk = k[2]-k[1]
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
