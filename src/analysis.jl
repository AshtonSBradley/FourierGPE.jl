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
