# test arrays

# array methods
linspace(a,b,n) = LinRange(a,b,n) |> collect

# function xvecs(Lx,Nx)
#     x = linspace(-Lx/2,Lx/2,Nx+1)[2:end]
#     return (x,)
# end
#
# function xvecs(Lx,Ly,Nx,Ny)
#     x = linspace(-Lx/2,Lx/2,Nx+1)[2:end]
#     y = linspace(-Ly/2,Ly/2,Ny+1)[2:end]
#     return x,y
# end
#
# function xvecs(Lx,Ly,Lz,Nx,Ny,Nz)
#     x = linspace(-Lx/2,Lx/2,Nx+1)[2:end]
#     y = linspace(-Ly/2,Ly/2,Ny+1)[2:end]
#     z = linspace(-Lz/2,Lz/2,Nz+1)[2:end]
#     return x,y,z
# end
#
# function kvecs(Lx,Nx)
#     @assert iseven(Nx)
#     nkx = 0:Int(Nx/2)
#     kx = [nkx[1:end-1];-reverse(nkx[2:end])]*2*π/Lx
#     return (kx,)
# end
#
# function kvecs(Lx,Ly,Nx,Ny)
#     @assert iseven(Nx)
#     @assert iseven(Ny)
#     nkx = 0:Int(Nx/2)
#     nky = 0:Int(Ny/2)
#     kx = [nkx[1:end-1];-reverse(nkx[2:end])]*2*π/Lx
#     ky = [nky[1:end-1];-reverse(nky[2:end])]*2*π/Ly
#     return kx, ky
# end
#
# function kvecs(Lx,Ly,Lz,Nx,Ny,Nz)
#     @assert iseven(Nx)
#     @assert iseven(Ny)
#     @assert iseven(Nz)
#     nkx = 0:Int(Nx/2); nky = 0:Int(Ny/2); nkz = 0:Int(Nz/2)
#     kx = [nkx[1:end-1];-reverse(nkx[2:end])]*2*π/Lx
#     ky = [nky[1:end-1];-reverse(nky[2:end])]*2*π/Ly
#     kz = [nkz[1:end-1];-reverse(nkz[2:end])]*2*π/Lz
#     return kx, ky, kz
# end





# function k2(Lx,Nx)
#     K = kvecs(Lx,Nx)
#     return [kx^2 for kx in K[1]] |> complex
# end
#
# function k2(Lx,Ly,Nx,Ny)
#     K = kvecs(Lx,Ly,Nx,Ny)
#     return [kx^2 + ky^2 for kx in K[1], ky in K[2]] |> complex
# end
#
# function k2(Lx,Ly,Lz,Nx,Ny,Nz)
#     K = kvecs(Lx,Ly,Lz,Nx,Ny,Nz)
#     return [kx^2 + ky^2 + kz^2 for kx in K[1], ky in K[2], kz in K[3]] |> complex
# end



# function makearrays(Lx,Nx)
#     X = xvecs(Lx,Nx)
#     K = kvecs(Lx,Nx)
#     x = X; kx = K
#     dX = diff(x[1])[1] |> Tuple
#     dK = diff(kx[1])[1] |> Tuple
#     return X,K,dX,dK
# end
#
# function makearrays(Lx,Ly,Nx,Ny)
#     X = xvecs(Lx,Ly,Nx,Ny)
#     K = kvecs(Lx,Ly,Nx,Ny)
#     x,y = X; kx,ky = K
#     dX = diff(x)[1],diff(y)[1]
#     dK = diff(kx[:])[1],diff(ky[:])[1]
#     return X,K,dX,dK
# end
#
# function makearrays(Lx,Ly,Lz,Nx,Ny,Nz)
#     X = xvecs(L,N)
#     K = kvecs(L,N)
#     x,y,z = X; kx,ky,kz = K
#     dX = diff(x)[1],diff(y)[1],diff(z)[1]
#     dK = diff(kx[:])[1],diff(ky[:])[1],diff(kz[:])[1]
#     return X,K,dX,dK
# end

# array methods
linspace(a,b,n) = LinRange(a,b,n) |> collect

function xvecs(L,N)
    X = []
    for (λ,ν) in zip(L,N)
        x = linspace(-λ/2,λ/2,ν+1)[2:end]
        push!(X,x)
    end
    return X |> Tuple
end

function kvecs(L,N)
    @assert all(iseven.(N))
    K = []
    for (λ,ν) in zip(L,N)
        nkx = 0:Int(ν/2)
        kx = [nkx[1:end-1];-reverse(nkx[2:end])]*2*π/λ
        push!(K,kx)
    end
    return K |> Tuple
end

function k2(L,N)
    M = length(L)
    K = kvecs(L,N)
    if M==1
    K2 = [k^2 for k ∈ K[1] ]
    elseif M==2
    K2 = [kx^2 + ky^2 for kx in K[1], ky in K[2]]
    elseif M==3
    K2 = [kx^2 + ky^2 + kz^2 for kx in K[1], ky in K[2], kz in K[3]]
    end
    return K2 |> complex
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
`A = crandn_array(M)`

Make `2x2` complex `randn()` array of dimension `M` (a placeholder)."""
function crandn_array(M)
    a = Int.(ones(M)).+1
    return randn(a...) |> complex
end



L = (10.,30.)
N = (30,20)

X = xvecsall(L,N)

K = kvecsall(L,N)

k3 = k2all(L,N)

k3test = k2(L...,N...)

all(k3test==k3)

L = (10.,30.,40.)
N = (30,20,10)
K = kvecs(L,N)

k3 = k2all(L,N)

k3test = k2(L...,N...)

all(k3test==k3)

L = (10.,)
N = (30,)
K = kvecs(L,N)

k3 = k2all(L,N)

k3test = k2(L...,N...)

all(k3test==k3)

K

function k2b(L,N)
    K = kvecs(L,N)
    k2 = zeros(N)
    for i in eachindex(N)
        for j in eachindex(K[i])
        k = K[i]
        k2
    for i ∈ eachindex(L)
        M = length(K[i])
        firstdims = []
        push!(firstdims,1)
        k = reshape(K[i][:],firstdims |> tuple ...,M)
    end
    return k2
end

reshape(k[:],length(k))

K2test =

k0 = rand(N...)
i = CartesianIndices(k0)
j = LinearIndices(k0)
k0[i[2,1,1]]
k0[j[2,1,1]]
for i ∈ CartesianIndices(k0)
    println(i)
end


L=(9.,8.)
N=(20,12)
X,K,dX,dK = makearrays(L...,N...)
Xtest,Ktest,dXtest,dKtest = makearraysall(L,N)

all(X==Xtest)
all(K==Ktest)
all(dK==dKtest)
all(dX==dXtest)


L=(9.,)
N=(20,)
X=xvecs(L,N)
x = X[1]
diff(x)[1]

X,K,dX,dK = makearrays(L...,N...)
Xtest,Ktest,dXtest,dKtest = makearraysall(L,N)

all(X==Xtest)
all(K==Ktest)
all(dK==dKtest)
all(dX==dXtest)


L=(9.,8.,12.)
N=(20,12,30)
X,K,dX,dK = makearrays(L...,N...)
Xtest,Ktest,dXtest,dKtest = makearraysall(L,N)

all(X==Xtest)
all(K==Ktest)
all(dK==dKtest)
all(dX==dXtest)


X = xvecs(L...,N...)
