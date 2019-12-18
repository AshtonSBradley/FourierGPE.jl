import Pkg;
Pkg.activate("FourierGPE")
ENV["JULIA_FFTW_PROVIDER"]="MKL"
Pkg.build("FFTW")
