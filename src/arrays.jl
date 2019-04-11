# array methods
linspace(a,b,n) = LinRange(a,b,n) |> collect

function xvecs(Lx,Nx)
    x = linspace(-Lx/2,Lx/2,Nx+1)[2:end]
    return (x,)
end

function kvecs(Lx,Nx)
    @assert iseven(Nx)
    nkx = 0:Int(Nx/2)
    kx = [nkx[1:end-1];-reverse(nkx[2:end])]*2*π/Lx
    return (kx,)
end

function xvecs(Lx,Ly,Nx,Ny)
    x = linspace(-Lx/2,Lx/2,Nx+1)[2:end]
    y = linspace(-Ly/2,Ly/2,Ny+1)[2:end]
    return x,y
end

function xvecs(Lx,Ly,Lz,Nx,Ny,Nz)
    x = linspace(-Lx/2,Lx/2,Nx+1)[2:end]
    y = linspace(-Ly/2,Ly/2,Ny+1)[2:end]
    z = linspace(-Lz/2,Lz/2,Nz+1)[2:end]
    return x,y,z
end

function kvecs(Lx,Ly,Nx,Ny)
    @assert iseven(Nx)
    @assert iseven(Ny)
    nkx = 0:Int(Nx/2)
    nky = 0:Int(Ny/2)
    kx = [nkx[1:end-1];-reverse(nkx[2:end])]*2*π/Lx
    ky = [nky[1:end-1];-reverse(nky[2:end])]*2*π/Ly
    return kx, ky
end

function kvecs(Lx,Ly,Lz,Nx,Ny,Nz)
    @assert iseven(Nx)
    @assert iseven(Ny)
    @assert iseven(Nz)
    nkx = 0:Int(Nx/2); nky = 0:Int(Ny/2); nkz = 0:Int(Nz/2)
    kx = [nkx[1:end-1];-reverse(nkx[2:end])]*2*π/Lx
    ky = [nky[1:end-1];-reverse(nky[2:end])]*2*π/Ly
    kz = [nkz[1:end-1];-reverse(nkz[2:end])]*2*π/Lz
    return kx, ky, kz
end

function k2(Lx,Nx)
    K = kvecs(Lx,Nx)
    return [kx^2 for kx in K[1]] |> complex
end

function k2(Lx,Ly,Nx,Ny)
    K = kvecs(Lx,Ly,Nx,Ny)
    return [kx^2 + ky^2 for kx in K[1], ky in K[2]] |> complex
end

function k2(Lx,Ly,Lz,Nx,Ny,Nz)
    K = kvecs(Lx,Ly,Lz,Nx,Ny,Nz)
    return [kx^2 + ky^2 + kz^2 for kx in K[1], ky in K[2], kz in K[3]] |> complex
end

function makearrays(Lx,Nx)
    X = (xvecs(Lx,Nx),)
    K = (kvecs(Lx,Nx),)
    x = X; kx = K
    dX = diff(x)[1]
    dK = diff(kx[:])[1]
    return X,K,dX,dK
end

function makearrays(Lx,Ly,Nx,Ny)
    X = xvecs(Lx,Ly,Nx,Ny)
    K = kvecs(Lx,Ly,Nx,Ny)
    x,y = X; kx,ky = K
    dX = diff(x)[1],diff(y)[1]
    dK = diff(kx[:])[1],diff(ky[:])[1]
    return X,K,dX,dK
end

function makearrays(Lx,Ly,Lz,Nx,Ny,Nz)
    X = xvecs(Lx,Ly,Lz,Nx,Ny,Lz)
    K = kvecs(Lx,Ly,Lz,Nx,Ny,Lz)
    x,y,z = X; kx,ky,kz = K
    dX = diff(x)[1],diff(y)[1],diff(z)[1]
    dK = diff(kx[:])[1],diff(ky[:])[1],diff(kz[:])[1]
    return X,K,dX,dK
end

"""
`A = crandn_array(M)`

Make `2x2` complex `randn()` array of dimension `M` (a placeholder)."""
function crandn_array(M)
    a = Int.(ones(M)).+1
    return randn(a...) |> complex
end
