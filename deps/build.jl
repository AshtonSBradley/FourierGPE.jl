import Pkg;
@info "Building FFTW#master with MKL"
ENV["JULIA_FFTW_PROVIDER"]="MKL"
Pkg.build("FFTW")
