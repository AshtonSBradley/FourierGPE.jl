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

function k2(L,N)
    M = length(L)
    K = kvecs(L,N)
    if M==1
        K2 = [k^2 for k ∈ K[1] ]
    elseif M==2
        K2 = [kx^2 + ky^2 for kx in K[1], ky in K[2]]
    elseif M==3
        K2 = [kx^2 + ky^2 + kz^2 for kx in K[1], ky in K[2], kz in K[3]]
    elseif M==4
        K2 = [kx^2 + ky^2 + kz^2 + kw^2 for kx in K[1], ky in K[2], kz in K[3], kw in K[4]]
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

# array methods
# linspace(a,b,n) = LinRange(a,b,n) |> collect
#
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
# function xvecsall(L,N)
#     X = []
#     for (λ,ν) in zip(L,N)
#         x = linspace(-λ/2,λ/2,ν+1)[2:end]
#         push!(X,x)
#     end
#     return X |> Tuple
# end
#
# function kvecsall(L,N)
#     @assert all(iseven.(N))
#     K = []
#     for (λ,ν) in zip(L,N)
#         nkx = 0:Int(ν/2)
#         kx = [nkx[1:end-1];-reverse(nkx[2:end])]*2*π/λ
#         push!(K,kx)
#     end
#     return X |> Tuple
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
#
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
#
# function makearrays(Lx,Nx)
#     X = (xvecs(Lx,Nx),)
#     K = (kvecs(Lx,Nx),)
#     x = X; kx = K
#     dX = diff(x)[1]
#     dK = diff(kx[:])[1]
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
#     X = xvecs(Lx,Ly,Lz,Nx,Ny,Lz)
#     K = kvecs(Lx,Ly,Lz,Nx,Ny,Lz)
#     x,y,z = X; kx,ky,kz = K
#     dX = diff(x)[1],diff(y)[1],diff(z)[1]
#     dK = diff(kx[:])[1],diff(ky[:])[1],diff(kz[:])[1]
#     return X,K,dX,dK
# end
#
# """
# `A = crandn_array(M)`
#
# Make `2x2` complex `randn()` array of dimension `M` (a placeholder)."""
# function crandn_array(M)
#     a = Int.(ones(M)).+1
#     return randn(a...) |> complex
# end
