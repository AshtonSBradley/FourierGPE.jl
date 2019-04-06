# transform methods

"""Measures that make fft, ifft 2-norm preserving."""
function dfft(x,k)
    dx = diff(x)[1]; dk = diff(k)[1]
    Dx = dx/sqrt(2*pi)
    Dk = length(k)*dk/sqrt(2*pi)
    return Dx, Dk
end

    function xspace(ψ,sim)
        @unpack T = sim
        return T.Tkx*ψ
    end

    function xspace!(ψ,sim)
        @unpack T = sim
        T.Tkx!*ψ
        return nothing
    end

    function kspace(ψ,sim)
        @unpack T = sim
        return T.Txk*ψ
    end

    function kspace!(ψ,sim)
        @unpack T = sim
        T.Txk!*ψ
        return nothing
    end

function maketransforms(Lx,Nx,Ly,Ny)
    x,y = xvecs(Lx,Nx,Ly,Ny)
    kx,ky,k2 = kvecs(Lx,Nx,Ly,Ny)

    #measures for Parseval tests
    Dx,Dkx = dfft(x,kx)
    Dy,Dky = dfft(y,ky)
    dx,dy  = diff(x)[1],diff(y)[1]
    dkx,dky = diff(kx)[1],diff(ky)[1]

    ψtest = one(x*y' |> complex)
    Txk = Dx*Dy*plan_fft(ψtest,flags=FFTW.MEASURE)
    ψtest = one(x*y' |> complex)
    Txk! = Dx*Dy*plan_fft!(ψtest,flags=FFTW.MEASURE)
    ψtest = one(x*y' |> complex)
    Tkx  = Dkx*Dky*plan_ifft(ψtest,flags=FFTW.MEASURE)
    ψtest = one(x*y' |> complex)
    Tkx!  = Dkx*Dky*plan_ifft!(ψtest,flags=FFTW.MEASURE)
    ψtest = one(x*y' |> complex)

return x,y,kx,ky,k2,dx,dy,dkx,dky,Dx,Dy,Dkx,Dky,Txk,Txk!,Tkx,Tkx!
end

@with_kw mutable struct Trans
    Txk::AbstractFFTs.ScaledPlan{Complex{Float64},FFTW.cFFTWPlan{Complex{Float64},-1,false,2},Float64} = 0.1*plan_fft(randn(2,2) |> complex)
    Txk!::AbstractFFTs.ScaledPlan{Complex{Float64},FFTW.cFFTWPlan{Complex{Float64},-1,true,2},Float64} = 0.1*plan_fft!(randn(2,2) |> complex)
    Tkx::AbstractFFTs.ScaledPlan{Complex{Float64},FFTW.cFFTWPlan{Complex{Float64},1,false,2},Float64} = 0.1*plan_ifft(randn(2,2) |> complex)
    Tkx!::AbstractFFTs.ScaledPlan{Complex{Float64},FFTW.cFFTWPlan{Complex{Float64},1,true,2},Float64} = 0.1*plan_ifft!(randn(2,2) |> complex)
end
