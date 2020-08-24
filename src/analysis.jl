"""
	gradient(psi::XField{D})

Compute the `D` vector gradient components of an `XField` of spatial dimension `D`.
The `D` gradient components returned are `D`-dimensional arrays.
"""
function gradient(psi::XField{1})
	@unpack psiX,K = psi; kx = K[1]; ψ = psiX
	ϕ = fft(ψ)
	ψx = ifft(im*kx.*ϕ)
	return ψx
end

function gradient(psi::XField{2})
	@unpack psiX,K = psi; kx,ky = K; ψ = psiX
	ϕ = fft(ψ)
	ψx = ifft(im*kx.*ϕ)
	ψy = ifft(im*ky'.*ϕ)
	return ψx,ψy
end

function gradient(psi::XField{3})
	@unpack psiX,K = psi; kx,ky,kz = K; ψ = psiX
	ϕ = fft(ψ)
	ψx = ifft(im*kx.*ϕ)
	ψy = ifft(im*ky'.*ϕ)
	ψz = ifft(im*reshape(kz,1,1,length(kz)).*ϕ)
	return ψx,ψy,ψz
end

"""
	current(psi::XField{D})

Compute the `D` current components of an `XField` of spatial dimension `D`.
The `D` cartesian components returned are `D`-dimensional arrays.
"""
function current(psi::XField{1})
	@unpack psiX,K = psi; kx = K[1]; ψ = psiX
	ψx = gradient(psi)
	jx = @. imag(conj(ψ)*ψx)
	return jx
end

function current(psi::XField{2})
	@unpack psiX,K = psi; kx,ky = K; ψ = psiX
    ψx,ψy = gradient(psi)
	jx = @. imag(conj(ψ)*ψx)
	jy = @. imag(conj(ψ)*ψy)
	return jx,jy
end

function current(psi::XField{3})
	@unpack psiX,K = psi; kx,ky,kz = K; ψ = psiX
    ψx,ψy,ψz = gradient(psi)
	jx = @. imag(conj(ψ)*ψx)
	jy = @. imag(conj(ψ)*ψy)
	jz = @. imag(conj(ψ)*ψz)
	return jx,jy,jz
end

"""
	velocity(psi::XField{D})

Compute the `D` velocity components of an `XField` of spatial dimension `D`.
The `D` velocities returned are `D`-dimensional arrays.
"""
function velocity(psi::XField{1})
	@unpack psiX,K = psi; kx = K[1]; ψ = psiX
	rho = abs2.(ψ)
    ψx = gradient(ψ)
	vx = @. imag(conj(ψ)*ψx)/rho; @. vx[isnan(vx)] = zero(vx[1])
	return vx
end

function velocity(psi::XField{2})
	@unpack psiX,K = psi; kx,ky = K; ψ = psiX
	rho = abs2.(ψ)
    ψx,ψy = gradient(psi)
	vx = @. imag(conj(ψ)*ψx)/rho; @. vx[isnan(vx)] = zero(vx[1])
	vy = @. imag(conj(ψ)*ψy)/rho; @. vy[isnan(vy)] = zero(vy[1])
	return vx,vy
end

function velocity(psi::XField{3})
	@unpack psiX,K = psi; kx,ky,kz = K; ψ = psiX
	rho = abs2.(ψ)
    ψx,ψy,ψz = gradient(psi)
	vx = @. imag(conj(ψ)*ψx)/rho; @. vx[isnan(vx)] = zero(vx[1])
	vy = @. imag(conj(ψ)*ψy)/rho; @. vy[isnan(vy)] = zero(vy[1])
	vz = @. imag(conj(ψ)*ψz)/rho; @. vz[isnan(vz)] = zero(vz[1])
	return vx,vy,vz
end

"""
	Wi,Wc = helmholtz(wx,wy,...,psi::XField{D})

Computes a 2 or 3 dimensional Helmholtz decomposition of the vector field with components
`wx`, `wy`, or `wx`, `wy`, `wz`. `psi` is passed to provide requisite arrays in `k`-space.
Returned fields `Wi`, `Wc` are incompressible and compressible respectively.
"""
function helmholtz(wx, wy, psi::XField{2})
    @unpack K, K2 = psi; kx, ky = K
    wxk = fft(wx); wyk = fft(wy)
    kdotw = @. kx * wxk + ky' * wyk
    wxkc = @. kdotw * kx / K2; wxkc[1] = zero(wxkc[1])
    wykc = @. kdotw * ky' / K2; wykc[1] = zero(wykc[1])
    wxki = @. wxk - wxkc
    wyki = @. wyk - wykc
    wxc = ifft(wxkc); wyc = ifft(wykc)
  	wxi = ifft(wxki); wyi = ifft(wyki)
  	Wi = (wxi, wyi); Wc = (wxc, wyc)
    return Wi, Wc
end

function helmholtz(wx, wy, wz, psi::XField{3})
    @unpack K, K2 = psi; kx, ky, kz = K
    wxk = fft(wx); wyk = fft(wy); wzk = fft(wz)
    kzt = reshape(kz, 1, 1, length(kz))
    kdotw = @. kx * wxk + ky' * wyk + kzt * wzk
    wxkc = @. kdotw * kx / K2; wxkc[1] = zero(wxkc[1])
    wykc = @. kdotw * ky' / K2; wykc[1] = zero(wykc[1])
    wzkc = @. kdotw * kzt / K2; wzkc[1] = zero(wzkc[1])
    wxki = @. wxk - wxkc
    wyki = @. wyk - wykc
    wzki = @. wzk - wzkc
    wxc = ifft(wxkc); wyc = ifft(wykc); wzc = ifft(wzkc)
    wxi = ifft(wxki); wyi = ifft(wyki); wzi = ifft(wzki)
  		Wi = (wxi, wyi, wzi); Wc = (wxc, wyc, wzc)
    return Wi, Wc
end

function helmholtz(W::NTuple{N,Array{Float64,N}}, psi::XField{N}) where N
    return helmholtz(W..., psi)
end

"""
	et,ei,ec = energydecomp(psi::Xfield{D})

Decomposes the hydrodynamic kinetic energy of `psi`, returning the total `et`, incompressible `ei`,
and compressible `ec` energy densities in position space. `D` can be 2 or 3 dimensions.
"""
function energydecomp(psi::XField{2})
    @unpack psiX = psi
    rho = abs2.(psiX)
    vx, vy = velocity(psi)
    wx = @. sqrt(rho)*vx; wy = @. sqrt(rho)*vy
    Wi, Wc = helmholtz(wx,wy,psi)
    wxi, wyi = Wi; wxc, wyc = Wc
    et = @. abs2(wx) + abs2(wy); et *= 0.5
    ei = @. abs2(wxi) + abs2(wyi); ei *= 0.5
    ec = @. abs2(wxc) + abs2(wyc); ec *= 0.5
    return et, ei, ec
end

function energydecomp(psi::XField{3})
	@unpack psiX = psi
    rho = abs2.(psiX)
    vx, vy, vz = velocity(psi)
    wx = @. sqrt(rho) * vx; wy = @. sqrt(rho) * vy; wz = @. sqrt(rho) * vz
    Wi, Wc = helmholtz(wx,wy,wz,psi)
    wxi, wyi, wzi = Wi; wxc, wyc, wzc = Wc
    et = @. abs2(wx) + abs2(wy) + abs2(wz); et *= 0.5
    ei = @. abs2(wxi) + abs2(wyi) + abs2(wzi); ei *= 0.5
    ec = @. abs2(wxc) + abs2(wyc) + abs2(wzc); ec *= 0.5
    return et, ei, ec
end

"""
	zeropad(A)

Zero-pad the 2D array `A` to twice the size with the same element type as `A`.
"""
function zeropad(a)
    s = size(a)
    (isodd.(s) |> any) && error("Array dims must be divisible by 2")
    S = @. 2 * s
    t = @. s / 2 |> Int
    z = Zeros{eltype(a)}(t...)
    M = Hcat([z; z], a, [z; z])
    return Vcat([z z z z], M, [z z z z]) |> Matrix
end

"""
	log10range(a,b,n)

Create a vector that is linearly spaced in log space, containing `n` values bracketed by `a` and `b`.
"""
function log10range(a,b,n)
	@assert a>0
    x = LinRange(log10(a),log10(b),n)
    return @. 10^x
end

@doc raw"""
	A = convolve(ψ1,ψ2,X,K)

Computes the convolution of two complex fields according to

```math
A(\rho) = \int d^2r\;\psi_1^*(r+\rho)\psi_2(r)
```
using FFTW.
"""
function convolve(ψ1,ψ2,X,K)
    n = length(X)
    DX,DK = dfftall(X,K)
	ϕ1 = zeropad(conj.(ψ1))
    ϕ2 = zeropad(ψ2)

	χ1 = fft(ϕ1)*prod(DX)
	χ2 = fft(ϕ2)*prod(DX)
	return ifft(χ1.*χ2)*prod(DK)*(2*pi)^(n/2) |> fftshift
end

@doc raw"""
	autocorrelate(ψ,X,K)

Return the autocorrelation integral of a complex field ``\psi``, ``A``, given by

```
A(\rho)=\int d^2r\;\psi^*(r-\rho)\psi(r)
```

defined on a cartesian grid on a cartesian grid using FFTW.

`X` and `K` are tuples of vectors `x`,`y`,`kx`, `ky`.

This method is useful for evaluating spectra from cartesian data.
"""
function autocorrelate(ψ,X,K)
    n = length(X)
    DX,DK = dfftall(X,K)
    ϕ = zeropad(ψ)
	χ = fft(ϕ)*prod(DX)
	return ifft(abs2.(χ))*prod(DK)*(2*pi)^(n/2) |> fftshift
end

function bessel_reduce(k,x,y,C)
    dx,dy = x[2]-x[1],y[2]-y[1]
    Nx = 2*length(x)
    Lx = x[end] - x[begin] + dx
    xp = LinRange(-Lx,Lx,Nx+1)[1:Nx]
    yp = xp
    ρ = hypot.(xp,yp')
    E = zero(k)
    @tullio E[i] = real(besselj0(k[i]*ρ[p,q])*C[p,q])
    @. E *= k*dx*dy/2/pi 
    return E 
end

function sinc_reduce(k,x,y,z,C)
    dx,dy,dz = x[2]-x[1],y[2]-y[1],z[2]-z[1]
    Nx = 2*length(x)    # assumes grids are same in each dimension
    Lx = x[end] - x[begin] + dx
    xp = LinRange(-Lx,Lx,Nx+1)[1:Nx]
    yp = xp'; zp = reshape(xp,1,1,Nx)
    ρ = hypot.(xp,yp,zp)
    E = zero(k)
    @tullio E[i] = real(sinc(k[i]*ρ[p,q,r]/π)*C[p,q,r]) # julia def. of sinc
    @. E *= k.^2*dx*dy*dz/2/pi^2  
    return E 
end

"""
	kinetic_edensity(k,ψ,X,K)

Calculates the kinetic enery spectrum for wavefunction ``\\psi``, at the
points `k`. Arrays `X`, `K` should be computed using `makearrays`.
"""
function kinetic_edensity(k,ψ::Array{Complex{Float64},2},X,K)
    x,y = X; kx,ky = K
    dx,dy = diff(x)[1],diff(y)[1]
    DX,DK = dfftall(X,K)
    k2field = k2(K)
    psi = XField(ψ,X,K,k2field)
    ψx,ψy = gradient(psi)

	cx = autocorrelate(ψx,X,K)
	cy = autocorrelate(ψy,X,K)
    C = @. 0.5(cx + cy)

    return bessel_reduce(k,x,y,C)
end

function kinetic_edensity(k,ψ::Array{Complex{Float64},3},X,K)
    x,y,z = X; kx,ky,kz = K
    dx,dy,dz = x[2]-x[1],y[2]-y[1],z[2]-z[1]
    DX,DK = dfftall(X,K)
    k2field = k2(K)
    psi = XField(ψ,X,K,k2field)
    ψx,ψy,ψz = gradient(psi)

	cx = autocorrelate(ψx,X,K)
    cy = autocorrelate(ψy,X,K)
    cz = autocorrelate(ψz,X,K)
    C = @. 0.5(cx + cy + cz)

    return sinc_reduce(k,x,y,z,C)
end

"""
	incompressible_cspectrum(k,ψ,X,K)

Caculate the incompressible velocity correlation spectrum for wavefunction ``\\psi``, via Helmholtz decomposition.
Input arrays `X`, `K` must be computed using `makearrays`.
"""
function incompressible_cspectrum(k,ψ::Array{Complex{Float64},2},X,K)
    x,y = X; kx,ky = K
 	dx,dy = x[2]-x[1],y[2]-y[1] 
	DX,DK = dfftall(X,K)
    k2field = k2(K)
    psi = XField(ψ,X,K,k2field)
    vx,vy = velocity(psi)
    wx = @. abs(ψ)*vx; wy = @. abs(ψ)*vy
    Wi, Wc = helmholtz(wx,wy,psi)
    wx,wy = Wi

	cx = autocorrelate(wx,X,K)
	cy = autocorrelate(wy,X,K)
    C = @. 0.5*(cx + cy)

    return bessel_reduce(k,x,y,C)
end

function incompressible_cspectrum(k,ψ::Array{Complex{Float64},3},X,K)
    x,y,z = X; kx,ky,kz = K
 	dx,dy,dz = x[2]-x[1],y[2]-y[1],z[2]-z[1]
	DX,DK = dfftall(X,K)
    k2field = k2(K)
    psi = XField(ψ,X,K,k2field)
    vx,vy,vz = velocity(psi)
    wx = @. abs(ψ)*vx; wy = @. abs(ψ)*vy; wz = @. abs(ψ)*vz
    Wi, Wc = helmholtz(wx,wy,wz,psi)
    wx,wy,wz = Wi

	cx = autocorrelate(wx,X,K)
    cy = autocorrelate(wy,X,K)
    cy = autocorrelate(wz,X,K)
    C = @. 0.5*(cx + cy + cz)

    return sinc_reduce(k,x,y,z,C)
end

"""
	compressible_cspectrum(k,ψ,X,K)

Caculate the compressible kinetic enery spectrum for wavefunction ``\\psi``, via Helmholtz decomposition.
Input arrays `X`, `K` must be computed using `makearrays`.
"""
function compressible_cspectrum(k,ψ::Array{Complex{Float64},2},X,K)
    x,y = X; kx,ky = K
 	dx,dy = x[2]-x[1],y[2]-y[1]
	DX,DK = dfftall(X,K)
    k2field = k2(K)
    psi = XField(ψ,X,K,k2field)
    vx,vy = velocity(psi)
    wx = @. abs(ψ)*vx; wy = @. abs(ψ)*vy
    Wi, Wc = helmholtz(wx,wy,psi)
    wx,wy = Wc

	cx = autocorrelate(wx,X,K)
	cy = autocorrelate(wy,X,K)
    C = @. 0.5*(cx + cy)

    return bessel_reduce(k,x,y,C)
end

function compressible_cspectrum(k,ψ::Array{Complex{Float64},3},X,K)
    x,y,z = X; kx,ky,kz = K
 	dx,dy,dz = x[2]-x[1],y[2]-y[1],z[2]-z[1]
	DX,DK = dfftall(X,K)
    k2field = k2(K)
    psi = XField(ψ,X,K,k2field)
    vx,vy,vz = velocity(psi)
    wx = @. abs(ψ)*vx; wy = @. abs(ψ)*vy; wz = @. abs(ψ)*vz
    Wi, Wc = helmholtz(wx,wy,wz,psi)
    wx,wy,wz = Wc

	cx = autocorrelate(wx,X,K)
    cy = autocorrelate(wy,X,K)
    cz = autocorrelate(wz,X,K)
    C = @. 0.5*(cx + cy + cz)

    return sinc_reduce(k,x,y,z,C)
end

"""
	qpressure_cspectrum(k,ψ,X,K)

Caculate the quantum pressure correlation spectrum for wavefunction ``\\psi``.
Input arrays `X`, `K` must be computed using `makearrays`.
"""
function qpressure_cspectrum(k,ψ::Array{Complex{Float64},2},X,K)
    x,y = X; kx,ky = K
    dx,dy = x[2]-x[1],y[2]-y[1]
	DX,DK = dfftall(X,K)
    k2field = k2(K)
    psi = XField(abs.(ψ) |> complex,X,K,k2field)
    wx,wy = gradient(psi)

	cx = autocorrelate(wx,X,K)
	cy = autocorrelate(wy,X,K)
    C = @. 0.5*(cx + cy)

    return bessel_reduce(k,x,y,C)
end

function qpressure_cspectrum(k,ψ::Array{Complex{Float64},3},X,K)
    x,y,z = X; kx,ky,kz = K
    dx,dy,dz = x[2]-x[1],y[2]-y[1],z[2]-z[1]
	DX,DK = dfftall(X,K)
    k2field = k2(K)
    psi = XField(abs.(ψ) |> complex,X,K,k2field)
    wx,wy,wz = gradient(psi)

	cx = autocorrelate(wx,X,K)
    cy = autocorrelate(wy,X,K)
    cz = autocorrelate(wz,X,K)
    C = @. 0.5*(cx + cy + cz)

    return sinc_reduce(k,x,y,z,C)
end

"""
    incompressible_edensity(k,ψ,X,K)

Calculates the kinetic energy density of the incompressible velocity field in the wavefunction ``\\psi``, at the
points `k`. Arrays `X`, `K` should be computed using `makearrays`.
"""
function incompressible_edensity(k,ψ::Array{Complex{Float64},2},X,K)
    x,y = X; kx,ky = K
    dx,dy = x[2]-x[1],y[2]-y[1] 
	DX,DK = dfftall(X,K)
    k2field = k2(K)
    psi = XField(ψ,X,K,k2field)
    vx,vy = velocity(psi)

    ux = @. abs(ψ)*vx; uy = @. abs(ψ)*vy 
    Wi, Wc = helmholtz(ux,uy,psi)
    wix,wiy = Wi
    U = @. exp(im*angle(ψ))
    @. wix *= U # restore phase factors
    @. wiy *= U

	cx = autocorrelate(wix,X,K)
	cy = autocorrelate(wiy,X,K)
    C = @. 0.5*(cx + cy)
    return bessel_reduce(k,x,y,C)
end

function incompressible_edensity(k,ψ::Array{Complex{Float64},3},X,K)
    x,y,z = X; kx,ky,kz = K
    dx,dy,dz = x[2]-x[1],y[2]-y[1],z[2]-z[1]
	DX,DK = dfftall(X,K)
    k2field = k2(K)
    psi = XField(ψ,X,K,k2field)
    vx,vy,vz = velocity(psi)

    ux = @. abs(ψ)*vx; uy = @. abs(ψ)*vy; uz = @. abs(ψ)*vz
    Wi, Wc = helmholtz(ux,uy,uz,psi)
    wix,wiy,wiz = Wi
    U = @. exp(im*angle(ψ))
    @. wix *= U # restore phase factors
    @. wiy *= U
    @. wiz *= U

	cx = autocorrelate(wix,X,K)
    cy = autocorrelate(wiy,X,K)
    cz = autocorrelate(wiz,X,K)
    C = @. 0.5*(cx + cy + cz)
    return sinc_reduce(k,x,y,z,C)
end

"""
    compressible_edensity(k,ψ,X,K)

Calculates the kinetic energy density of the compressible velocity field in the wavefunction ``\\psi``, at the
points `k`. Arrays `X`, `K` should be computed using `makearrays`.
"""
function compressible_edensity(k,ψ::Array{Complex{Float64},2},X,K)
    x,y = X; kx,ky = K
    dx,dy = x[2]-x[1],y[2]-y[1] 
	DX,DK = dfftall(X,K)
    k2field = k2(K)
    psi = XField(ψ,X,K,k2field)
    vx,vy = velocity(psi)

    ux = @. abs(ψ)*vx; uy = @. abs(ψ)*vy 
    Wi, Wc = helmholtz(ux,uy,psi)
    wcx,wcy = Wc
    U = @. exp(im*angle(ψ))
    @. wcx *= U # restore phase factors
    @. wcy *= U

	cx = autocorrelate(wcx,X,K)
	cy = autocorrelate(wcy,X,K)
    C = @. 0.5*(cx + cy)
    return bessel_reduce(k,x,y,C)
end

function compressible_edensity(k,ψ::Array{Complex{Float64},3},X,K)
    x,y,z = X; kx,ky,kz = K
    dx,dy,dz = x[2]-x[1],y[2]-y[1],z[2]-z[1]
	DX,DK = dfftall(X,K)
    k2field = k2(K)
    psi = XField(ψ,X,K,k2field)
    vx,vy,vz = velocity(psi)

    ux = @. abs(ψ)*vx; uy = @. abs(ψ)*vy; uz = @. abs(ψ)*vz
    Wi, Wc = helmholtz(ux,uy,uz,psi)
    wcx,wcy,wcz = Wc
    U = @. exp(im*angle(ψ))
    @. wcx *= U # restore phase factors
    @. wcy *= U
    @. wcz *= U

	cx = autocorrelate(wcx,X,K)
    cy = autocorrelate(wcy,X,K)
    cz = autocorrelate(wcz,X,K)
    C = @. 0.5*(cx + cy + cz)
    return sinc_reduce(k,x,y,z,C)
end

"""
    qpressure_edensity(k,ψ,X,K)

Energy density of the quantum pressure in the wavefunction ``\\psi``, at the
points `k`. Arrays `X`, `K` should be computed using `makearrays`.
"""
function qpressure_edensity(k,ψ::Array{Complex{Float64},2},X,K)
    x,y = X; kx,ky = K
    dx,dy = x[2]-x[1],y[2]-y[1] 
	DX,DK = dfftall(X,K)
    k2field = k2(K)
    psi = XField(abs.(ψ) |> complex,X,K,k2field)
    rnx,rny = gradient(psi)
    U = @. exp(im*angle(ψ))
    @. rnx *= U # restore phase factors
    @. rny *= U 

	cx = autocorrelate(rnx,X,K)
	cy = autocorrelate(rny,X,K)
    C = @. 0.5*(cx + cy)
    return bessel_reduce(k,x,y,C)
end

function qpressure_edensity(k,ψ::Array{Complex{Float64},3},X,K)
    x,y,z = X; kx,ky,kz = K
    dx,dy,dz = x[2]-x[1],y[2]-y[1],z[2]-z[1]
	DX,DK = dfftall(X,K)
    k2field = k2(K)
    psi = XField(abs.(ψ) |> complex,X,K,k2field)
    rnx,rny,rnz = gradient(psi)
    U = @. exp(im*angle(ψ))
    @. rnx *= U # restore phase factors
    @. rny *= U 
    @. rnz *= U 

	cx = autocorrelate(rnx,X,K)
    cy = autocorrelate(rny,X,K)
    cz = autocorrelate(rnz,X,K)
    C = @. 0.5*(cx + cy + cz)
    return sinc_reduce(k,x,y,z,C)
end

## coupling terms

"""
    ic_edensity(k,ψ,X,K)

Energy density of the incompressible-compressible interaction in the wavefunction ``\\psi``, at the
points `k`. Arrays `X`, `K` should be computed using `makearrays`.
"""
function ic_edensity(k,ψ::Array{Complex{Float64},2},X,K)
    x,y = X; kx,ky = K
    dx,dy = x[2]-x[1],y[2]-y[1] 
	DX,DK = dfftall(X,K)
    k2field = k2(K)
    psi = XField(ψ,X,K,k2field)
    vx,vy = velocity(psi)

    ux = @. abs(ψ)*vx; uy = @. abs(ψ)*vy 
    Wi, Wc = helmholtz(ux,uy,psi)
    wix,wiy = Wi; wcx,wcy = Wc
    U = @. exp(im*angle(ψ))
    @. wix *= im*U # restore phase factors and make u -> w fields
    @. wiy *= im*U
    @. wcx *= im*U 
    @. wcy *= im*U

    cicx = convolve(wix,wcx,X,K) 
    ccix = convolve(wcx,wix,X,K)
    cicy = convolve(wiy,wcy,X,K) 
    cciy = convolve(wcy,wiy,X,K)
    C = @. 0.5*(cicx + ccix + cicy + cciy)  
    return bessel_reduce(k,x,y,C)
end

function ic_edensity(k,ψ::Array{Complex{Float64},3},X,K)
    x,y,z = X; kx,ky,kz = K
    dx,dy,dz = x[2]-x[1],y[2]-y[1],z[2]-z[1]
	DX,DK = dfftall(X,K)
    k2field = k2(K)
    psi = XField(ψ,X,K,k2field)
    vx,vy,vz = velocity(psi)

    ux = @. abs(ψ)*vx; uy = @. abs(ψ)*vy; uz = @. abs(ψ)*vz 
    Wi, Wc = helmholtz(ux,uy,uz,psi)
    wix,wiy,wiz = Wi; wcx,wcy,wcz = Wc
    U = @. exp(im*angle(ψ))
    @. wix *= im*U # restore phase factors and make u -> w fields
    @. wiy *= im*U
    @. wiz *= im*U   
    @. wcx *= im*U 
    @. wcy *= im*U
    @. wcz *= im*U

    cicx = convolve(wix,wcx,X,K) 
    ccix = convolve(wcx,wix,X,K)
    cicy = convolve(wiy,wcy,X,K) 
    cciy = convolve(wcy,wiy,X,K)
    cicz = convolve(wiz,wcz,X,K) 
    cciz = convolve(wcz,wiz,X,K)
    C = @. 0.5*(cicx + ccix + cicy + cciy + cicz + cciz)  
    return sinc_reduce(k,x,y,z,C)
end

"""
    iq_edensity(k,ψ,X,K)

Energy density of the incompressible-quantum pressure interaction in the wavefunction ``\\psi``, at the
points `k`. Arrays `X`, `K` should be computed using `makearrays`.
"""
function iq_edensity(k,ψ::Array{Complex{Float64},2},X,K)
    x,y = X; kx,ky = K
    dx,dy = x[2]-x[1],y[2]-y[1] 
	DX,DK = dfftall(X,K)
    k2field = k2(K)
    psi = XField(ψ,X,K,k2field)
    vx,vy = velocity(psi)

    ux = @. abs(ψ)*vx; uy = @. abs(ψ)*vy 
    Wi, Wc = helmholtz(ux,uy,psi)
    wix,wiy = Wi 

    psi = XField(abs.(ψ) |> complex,X,K,k2field)
    wqx,wqy = gradient(psi)

    U = @. exp(im*angle(ψ))
    @. wix *= im*U # phase factors and make u -> w fields
    @. wiy *= im*U
    @. wqx *= U
    @. wqy *= U

    ciqx = convolve(wix,wqx,X,K) 
    cqix = convolve(wqx,wix,X,K) 
    ciqy = convolve(wiy,wqy,X,K) 
    cqiy = convolve(wqy,wiy,X,K) 
    C = @. 0.5*(ciqx + cqix + ciqy + cqiy) 
    return bessel_reduce(k,x,y,C)
end

function iq_edensity(k,ψ::Array{Complex{Float64},3},X,K)
    x,y,z = X; kx,ky,kz = K
    dx,dy,dz = x[2]-x[1],y[2]-y[1],z[2]-z[1]
	DX,DK = dfftall(X,K)
    k2field = k2(K)
    psi = XField(ψ,X,K,k2field)
    vx,vy,vz = velocity(psi)

    ux = @. abs(ψ)*vx; uy = @. abs(ψ)*vy; uz = @. abs(ψ)*vz
    Wi, Wc = helmholtz(ux,uy,uz,psi)
    wix,wiy,wiz = Wi; 

    psi = XField(abs.(ψ) |> complex,X,K,k2field)
    wqx,wqy,wqz = gradient(psi)

    U = @. exp(im*angle(ψ))
    @. wix *= im*U # phase factors and make u -> w fields
    @. wiy *= im*U
    @. wiz *= im*U
    @. wqx *= U
    @. wqy *= U
    @. wqz *= U

    ciqx = convolve(wix,wqx,X,K) 
    cqix = convolve(wqx,wix,X,K) 
    ciqy = convolve(wiy,wqy,X,K) 
    cqiy = convolve(wqy,wiy,X,K) 
    ciqz = convolve(wiz,wqz,X,K) 
    cqiz = convolve(wqz,wiz,X,K) 
    C = @. 0.5*(ciqx + cqix + ciqy + cqiy + ciqz + cqiz) 
    return sinc_reduce(k,x,y,z,C)
end


"""
    cq_edensity(k,ψ,X,K)

Energy density of the compressible-quantum pressure interaction in the wavefunction ``\\psi``, at the
points `k`. Arrays `X`, `K` should be computed using `makearrays`.
"""
function cq_edensity(k,ψ::Array{Complex{Float64},2},X,K)
    x,y = X; kx,ky = K
    dx,dy = x[2]-x[1],y[2]-y[1] 
	DX,DK = dfftall(X,K)
    k2field = k2(K)
    psi = XField(ψ,X,K,k2field)
    vx,vy = velocity(psi)

    ux = @. abs(ψ)*vx; uy = @. abs(ψ)*vy 
    Wi, Wc = helmholtz(ux,uy,psi)
    wcx,wcy = Wc 

    psi = XField(abs.(ψ) |> complex,X,K,k2field)
    wqx,wqy = gradient(psi)

    U = @. exp(im*angle(ψ))
    @. wcx *= im*U # phase factors and make u -> w fields
    @. wcy *= im*U
    @. wqx *= U
    @. wqy *= U

    ccqx = convolve(wcx,wqx,X,K) 
    cqcx = convolve(wqx,wcx,X,K) 
    ccqy = convolve(wcy,wqy,X,K) 
    cqcy = convolve(wqy,wcy,X,K) 
    C = @. 0.5*(ccqx + cqcx + ccqy + cqcy) 
    return bessel_reduce(k,x,y,C)
end

function cq_edensity(k,ψ::Array{Complex{Float64},3},X,K)
    x,y,z = X; kx,ky,kz = K
    dx,dy,dz = x[2]-x[1],y[2]-y[1],z[2]-z[1]
	DX,DK = dfftall(X,K)
    k2field = k2(K)
    psi = XField(ψ,X,K,k2field)
    vx,vy,vz = velocity(psi)

    ux = @. abs(ψ)*vx; uy = @. abs(ψ)*vy; uz = @. abs(ψ)*vz
    Wi, Wc = helmholtz(ux,uy,uz,psi)
    wcx,wcy,wcz = Wc  

    psi = XField(abs.(ψ) |> complex,X,K,k2field)
    wqx,wqy,wqz = gradient(psi)

    U = @. exp(im*angle(ψ))
    @. wcx *= im*U # phase factors and make u -> w fields
    @. wcy *= im*U
    @. wcz *= im*U
    @. wqx *= U
    @. wqy *= U
    @. wqz *= U

    ccqx = convolve(wcx,wqx,X,K) 
    cqcx = convolve(wqx,wcx,X,K) 
    ccqy = convolve(wcy,wqy,X,K) 
    cqcy = convolve(wqy,wcy,X,K) 
    ccqz = convolve(wcz,wqz,X,K) 
    cqcz = convolve(wqz,wcz,X,K) 
    C = @. 0.5*(ccqx + cqcx + ccqy + cqcy + ccqz + cqcz) 
    return sinc_reduce(k,x,y,z,C)
end