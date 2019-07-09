<<<<<<< HEAD
# setup types for analysis
using Parameters, FFTW, AbstractFFTs

abstract type Field end

struct XField{D} <: Field
    PsiX::Array{Complex{Float64},D}
    X::NTuple{D}
    K::NTuple{D}
    K2::Array{Float64,D}
end

struct KField{D} <: Field
    PsiK::Array{Complex{Float64},D}
    X::NTuple{D}
    K::NTuple{D}
    K2::Array{Float64,D}
end

xvec(L,N) = LinRange(-L/2,L/2,N+1)[2:end] |> collect

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

import AbstractFFTs.fft, AbstractFFTs.ifft
fft(psi::XField) = KField(fft(psi.PsiX),psi.X,psi.K,psi.K2)
ifft(psi::KField) = XField(ifft(psi.PsiK),psi.X,psi.K,psi.K2)

N = 100
psi = randn(N,N) + im*randn(N,N)

X = xvecs((1,2),(N,N))
K = kvecs((1,2),(N,N))
K2 = k2(K)

typeof(K)

xpsi = XField(psi,X,K,K2)
kpsi = KField(fft(psi),X,K,K2)
kpsi2 = fft(xpsi)


function velocity(psi::XField{2})
    @unpack PsiX, X,K = psi; kx,ky = K
    rho = abs2.(PsiX)
    ϕ = fft(PsiX)
    ψx = ifft(im*kx.*ϕ)
    ψy = ifft(im*ky'.*ϕ)
    vx = @. real(conj(PsiX)*ψx)/rho
    vy = @. real(conj(PsiX)*ψy)/rho
    return vx,vy
end

vx,vy = velocity(xpsi)

function velocity(phi::KField{2})
    @unpack KPsi, KGrids, k2 = phi; kx,ky = KGrids
    rho = abs2.(ψ)
    ϕ = fft(ψ)
    ψx = ifft(im*kx.*ϕ)
    ψy = ifft(im*ky.*ϕ)
    vx = @. real(conj(ψ)*ψx))/rho
    vy = @. real(conj(ψ)*ψy))/rho
    return vx,vy
end

function velocity(ψ,kx,ky,kz,k2)
    rho = abs2.(ψ)
    ϕ = fft(ψ)
    ψx = ifft(im*kx.*ϕ)
    ψy = ifft(im*ky.*ϕ)
    ψz = ifft(im*kz.*ϕ)
    vx = @. real(conj(ψ)*ψx))/rho
    vy = @. real(conj(ψ)*ψy))/rho
    vz = @. real(conj(ψ)*ψz))/rho
    return vx,vy,vz
end

function vorticity(ψ,kx,ky,k2)
    vx,vy = velocity(ψ,kx,ky,k2)
    vxk = fft(vx)
    vxy = ifft(im*ky.*vxk)
    vyk = fft(vy)
    vyx = ifft(im*ky.*vyk)
    w = @. vxy - vyx
    return w
end

function vorticity(ψ,kx,ky,kz,k2)
    vx,vy,vz = velocity(ψ,kx,ky,kz,k2)
    #TODO
    return wx,wy,wz
end

function energydecomp(ψ,kx,ky,kz,k2)
    rho = abs2.(ψ)
    vx,vy,vz = velocity(ψ,kx,ky,kz,k2)
    wx = @. sqrt(rho)*vx; wy = @. sqrt(rho)*vy; wz = @. sqrt(rho)*vx
    wxi,wyi,wzi,wxc,wyc,wzc = helmholtzdecomp(wx,wy,wz,kx,ky,kz,k2)
    et = @. abs2(wx)+abs2(wy)+abs2(wz); et *= 0.5
    ei = @. abs2(wxi)+abs2(wyi)+abs2(wzi); ei *= 0.5
    ec = @. abs2(wxc)+abs2(wyc)+abs2(wzc); ec *= 0.5
    return et,ei,ec
end

function energydecomp(ψ,kx,ky,k2)
    rho = abs2.(ψ)
    vx, vy = velocity(ψ,kx,ky,k2)
    wx = @. sqrt(rho)*vx; wy = @. sqrt(rho)*vy
    wxi, wyi, wxc, wyc = helmholtzdecomp(wx,wy,kx,ky,k2)
    et = @. abs2(wx)+abs2(wy); et *= 0.5
    ei = @. abs2(wxi)+abs2(wyi); ei *= 0.5
    ec = @. abs2(wxc)+abs2(wyc); ec *= 0.5
    return et,ei,ec
end

function helmholtzdecomp(wx,wy,kx,ky,k2)
    wxk = fft(wx); wyk = fft(wy)
    kdotw = @. kx*wxk + ky*wyk
    wxkc = @. kdotw*kx/k2; wxkc[1,1] = zero(wxkc[1,1]) #TODO necessary?
    wykc = @. kdotw*ky/k2; wykc[1,1] = zero(wykc[1,1])
    wxki = @. wxk - wxkc
    wyki = @. wyk - wykc
    wxc = ifft(wxkc); wyc = ifft(wykc)
    wxi = ifft(wxki); wyi = ifft(wyki)
    return wxi, wyi, wxc, wyc
end

function helmholtzdecomp(wx,wy,wz,kx,ky,kz,k2)
    #TODO check k sizes and do some careful tests!
    # ky = ky'
    # kz = reshape(kz,1,1,length(ky))
    wxk = fft(wx); wyk = fft(wy); wzk = fft(wz);
    kdotw = @. kx*wxk + ky*wyk + kz*wzk
    wxkc = @. kdotw*kx/k2; wxkc[1,1] = zero(wxkc[1,1])
    wykc = @. kdotw*ky/k2; wykc[1,1] = zero(wykc[1,1])
    wzkc = @. kdotw*kz/k2; wzkc[1,1] = zero(wzkc[1,1])
    wxki = @. wxk - wxkc
    wyki = @. wyk - wykc
    wzki = @. wzk - wzkc
    wxc = ifft(wxkc); wyc = ifft(wykc); wzc = ifft(wzkc)
    wxi = ifft(wxki); wyi = ifft(wyki); wzi = ifft(wzki)
    return wxi, wyi, wzi, wxc, wyc, wzc
end
=======
function velocity(psi::XField{2})
	@unpack psiX,K,K2 = psi; kx,ky=K; ψ = psiX
	rho = abs2.(ψ)
	ϕ = fft(ψ)
	ψx = ifft(im*kx.*ϕ)
	ψy = ifft(im*ky'.*ϕ)
	vx = @. imag(conj(ψ)*ψx)/rho
	vy = @. imag(conj(ψ)*ψy)/rho
	return vx,vy
end

function velocity(psi::XField{3})
	@unpack psiX,K,K2 = psi;kx,ky,kz=K; ψ = psiX
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

	function energydecomp(psi::XField{2})
		@unpack psiX,K,K2 = psi; kx,ky = K; ψ = psiX
	    rho = abs2.(ψ)
	    vx,vy = velocity(psi)
	    wx = @. sqrt(rho)*vx; wy = @. sqrt(rho)*vy
	    Wi,Wc = helmholtzdecomp(W,K,k2)
		wxi,wyi = Wi; wxc,wyc = Wc
	    et = @. abs2(wx)+abs2(wy); et *= 0.5
	    ei = @. abs2(wxi)+abs2(wyi); ei *= 0.5
	    ec = @. abs2(wxc)+abs2(wyc); ec *= 0.5
	    return et,ei,ec
	end

	function energydecomp(psi::XField{3})
		@unpack psiX,K,K2 = psi; kx,ky,kz = K; ψ = psiX
		    rho = abs2.(ψ)
		    vx,vy,vz = velocity(psi)
		    wx = @. sqrt(rho)*vx; wy = @. sqrt(rho)*vy; wz = @. sqrt(rho)*vz
		    Wi,Wc = helmholtzdecomp(W,K,k2)
			wxi,wyi,wzi = Wi; wxc,wyc,wzc = Wc
		    et = @. abs2(wx)+abs2(wy)+abs2(wz); et *= 0.5
		    ei = @. abs2(wxi)+abs2(wyi)+abs2(wzi); ei *= 0.5
		    ec = @. abs2(wxc)+abs2(wyc)+abs2(wzc); ec *= 0.5
		    return et,ei,ec
		end

	function helmholtzdecomp(W::NTuple{2},K::NTuple{2},K2)
		wx,wy = W; kx,ky = K
	    wxk = fft(wx); wyk = fft(wy)
	    kdotw = @. kx*wxk + ky'*wyk
	    wxkc = @. kdotw*kx/K2; wxkc[1,1] = zero(wxkc[1,1])
	    wykc = @. kdotw*ky/K2; wykc[1,1] = zero(wykc[1,1])
	    wxki = @. wxk - wxkc
	    wyki = @. wyk - wykc
	    wxc = ifft(wxkc); wyc = ifft(wykc)
	    wxi = ifft(wxki); wyi = ifft(wyki)
		Wi = (wxi,wyi); Wc = (wxc,wyc)
	    return Wi,Wc
	end

	function helmholtzdecomp(W::NTuple{3},K::NTuple{3},K2)
		wx,wy,wz = W; kx,ky,kz = K
	    wxk = fft(wx); wyk = fft(wy); wzk = fft(wz);
	    kdotw = @. kx*wxk + ky'*wyk + reshape(kz,1,1,length(kz))*wzk
	    wxkc = @. kdotw*kx/K2; wxkc[1,1] = zero(wxkc[1,1])
	    wykc = @. kdotw*ky/K2; wykc[1,1] = zero(wykc[1,1])
	    wzkc = @. kdotw*kz/K2; wzkc[1,1] = zero(wzkc[1,1])
	    wxki = @. wxk - wxkc
	    wyki = @. wyk - wykc
	    wzki = @. wzk - wzkc
	    wxc = ifft(wxkc); wyc = ifft(wykc); wzc = ifft(wzkc)
	    wxi = ifft(wxki); wyi = ifft(wyki); wzi = ifft(wzki)
		Wi = (wxi,wyi,wzi); Wc = (wxc,wyc,wzc)
	    return Wi,Wc
	end
>>>>>>> 31836ae5986f7dd75d608ba60861281fb2505075
