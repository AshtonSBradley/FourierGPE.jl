# array methods
linspace(a,b,n) = LinRange(a,b,n) |> collect

function xvecs(Lx,Ly,Nx,Ny)
    x = linspace(-Lx/2,Lx/2,Nx+1)[2:end]
    y = linspace(-Ly/2,Ly/2,Ny+1)[2:end]
    return x,y
end

function makearrays(Lx,Ly,Nx,Ny)
    x,y = xvecs(Lx,Ly,Nx,Ny)
    kx,ky,k2 = kvecs(Lx,Ly,Nx,Ny)
    dx,dy = diff(x)[1],diff(y)[1]
    dkx,dky = diff(kx[:])[1],diff(ky[:])[1]
    return x,y,kx,ky,dx,dy,dkx,dky
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
    k2 = [kx^2 + ky^2 for kx in kx, ky in ky] |> complex
    return kx, ky, k2
end

function kvecs(Lx,Ly,Lz,Nx,Ny,Nz)
    @assert iseven(Nx)
    @assert iseven(Ny)
    @assert iseven(Nz)
    nkx = 0:Int(Nx/2); nky = 0:Int(Ny/2); nkz = 0:Int(Nz/2)
    kx = [nkx[1:end-1];-reverse(nkx[2:end])]*2*π/Lx
    ky = [nky[1:end-1];-reverse(nky[2:end])]*2*π/Ly
    kz = [nkz[1:end-1];-reverse(nkz[2:end])]*2*π/Lz
    k2 = [kx^2 + ky^2 + kz^2 for kx in kx, ky in ky, kz in kz] |> complex
    return kx, ky, kz, k2
end
