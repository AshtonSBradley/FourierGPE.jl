"""
	vx,... = velocity(psi::XField{D})

Retruns the `D` velocity components of an `XField` of spatial dimension `D`. Velocities are `D`-dimensional arrays.
"""
function velocity(psi::XField{1})
	@unpack psiX,K = psi; kx = K[1]; ψ = psiX
	rho = abs2.(ψ)
	ϕ = fft(ψ)
	ψx = ifft(im*kx.*ϕ)
	vx = @. imag(conj(ψ)*ψx)/rho
	return vx
end

function velocity(psi::XField{2})
	@unpack psiX,K = psi; kx,ky = K; ψ = psiX
	rho = abs2.(ψ)
	ϕ = fft(ψ)
	ψx = ifft(im*kx.*ϕ)
	ψy = ifft(im*ky'.*ϕ)
	vx = @. imag(conj(ψ)*ψx)/rho
	vy = @. imag(conj(ψ)*ψy)/rho
	return vx,vy
end

function velocity(psi::XField{3})
	@unpack psiX,K = psi; kx,ky,kz = K; ψ = psiX
	rho = abs2.(ψ)
	ϕ = fft(ψ)
	ψx = ifft(im*kx.*ϕ)
	ψy = ifft(im*ky'.*ϕ)
	ψz = ifft(im*reshape(kz,1,1,length(kz)).*ϕ)
	vx = @. imag(conj(ψ)*ψx)/rho
	vy = @. imag(conj(ψ)*ψy)/rho
	vz = @. imag(conj(ψ)*ψz)/rho
	return vx,vy,vz
end

"""
	Wi,Wc = helmholtz(wx,wy,...,psi::XField{D})

Computes a 2 or 3 dimensional Helmholtz decomposition of the vector field with components
`wx`, `wy`, or `wx`, `wy`, `wz`. `psi` is passed to provide requisite arraysin `k` space.
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

Decomposes the kinetic energy of `psi`, returning the total `et`, incompressible `ei`,
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

@doc raw"""
	A = zeropad(a)

Zero-pad the 2D array `a` to twice the size with the same element type as `a`.
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

@doc raw"""
	x = log10range(a,b,n)

Create a logarithmically spaced vector. `x` is linear in log space, containing `n` values bracketed by `a` and `b`.
"""
function log10range(a,b,n)
    x = LinRange(log10(a),log10(b),n)
    return @. 10^x
end

@doc raw"""
	A = convolve(ψ1,ψ2,X,K)

Computes the convolution of two complex fields according to

``
A(\rho) = \int d^2r\;\psi_1(r-\rho)\psi_2(r)
``

using FFTW.
"""
function convolve(ψ1,ψ2,X,K)
	DX,DK = dfftall(X,K)
	ϕ1 = zeropad(ψ1)
	ϕ2 = zeropad(ψ2)
	χ1 = fft(ϕ1)*prod(DX)
	χ2 = fft(ϕ2)*prod(DX)
	return ifft(χ1.*χ2)*prod(DK) |> fftshift
end

@doc raw"""
	A = autocorrelate(ψ,X,K)

Evaluates the autocorrelation integral of a complex field ``\psi``.

``
A(\rho)=\int d^2r\;\psi^*(r-\rho)\psi(r)
``

defined on a cartesian grid on a cartesian grid using FFTW. `X` and `K` are tuples of vectors for `x`,`y`,`kx`, `ky`.

This method is useful for evaluating spectra from cartesian data.
"""
function autocorrelate(ψ,X,K)
	DX,DK = dfftall(X,K)
	ϕ = zeropad(ψ)
	χ = fft(ϕ)*prod(DX)
	return ifft(abs2.(χ))*prod(DK) |> fftshift
end

@doc raw"""
	Ek = kespectrum(k,ψ,X,K)

Caculates the kinetic enery spectrum for wavefunction `\psi`.
"""
function kespectrum(k,ψ,X,K)
	x,y = X; kx,ky = K
    dx,dy = diff(x)[1],diff(y)[1]
	DX,DK = dfftall(X,K)
    Nx = 2*length(x)
    Lx = x[end]-x[1] + dx
    xp = LinRange(-Lx,Lx,Nx+1)[1:Nx]
    yp = xp
    ρ = @. sqrt(xp^2 + yp'^2)
    # ψc = zeropad(ψ)
    # ϕc = fft(ψc)*prod(DX)
    # A = ifft(abs2.(ϕc))*prod(DK) |> fftshift
	A = autocorrelate(ψ,X,K)

    Ek = zero(k)
    for i in eachindex(k)
        κ = k[i]
        Ek[i]  = 0.5*κ^3*sum(@. besselj0(κ*ρ)*A)*dx*dy |> real
    end
    return Ek
end

@doc raw"""
	Ek = ikespectrum(k,ψ,X,K)

Caculates the incompressible kinetic enery spectrum for wavefunction `\psi`, using Helmholtz decomposition.
"""
function ikespectrum(k,ψ,X,K)
    x,y = X; kx,ky = K
 	dx,dy = diff(x)[1],diff(y)[1]
	DX,DK = dfftall(X,K)
    k2field = k2(K)
    psi = XField(ψ,X,K,k2field)
    vx,vy = velocity(psi)
    rho = abs2.(ψ)
    wx = @. sqrt(rho)*vx; wy = @. sqrt(rho)*vy
    Wi, Wc = helmholtz(wx,wy,psi)
    wix,wiy = Wi

    # transforms
    # Wix = zeropad(wix)
    # Wiy = zeropad(wiy)
    # Wixk = fft(Wix)*prod(DX)
    # Wiyk = fft(Wiy)*prod(DX)
    # @test sum(abs2.(Wix))*dx*dy ≈ sum(abs2.(Wixk))*(dkx/2)*(dky/2)

    # convolutions
    # Cix = ifft(abs2.(Wixk))*prod(DK) |> fftshift
    # Ciy = ifft(abs2.(Wiyk))*prod(DK) |> fftshift
	cwix = autocorrelate(wix,X,K)
	cwiy = autocorrelate(wiy,X,K)
    Ci = cwix .+ cwiy

    Nx = 2*length(x)
    Lx = x[end]-x[1] + dx
    xp = LinRange(-Lx,Lx,Nx+1)[1:Nx]
    yp = xp
    ρ = @. sqrt(xp^2 + yp'^2)

    Eki = zero(k)
    for i in eachindex(k)
        κ = k[i]
        Eki[i]  = 0.5*κ*sum(@. besselj0(κ*ρ)*Ci)*dx*dy |> real
    end
    return Eki
end
