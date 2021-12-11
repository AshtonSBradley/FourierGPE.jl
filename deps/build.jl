# MKL now supports 3D transforms, so...
# ENV["JULIA_FFTW_PROVIDER"]="FFTW"
# ENV["JULIA_FFTW_PROVIDER"]="MKL"
@info "Setting fft provider to MKL..."
FFTW.set_provider!("mkl")

