using Plots, LaTeXStrings, Pkg, Revise
gr(titlefontsize=12,size=(500,300),transpose=true,colorbar=false)

using FourierGPE
FFTW.set_num_threads(Sys.CPU_THREADS)

# ==== Units: ========================
# this example works in oscillator units
# convenient plot
function showpsi(ψ,x,y)
    p1 = heatmap(x,y,abs2.(ψ),aspectratio=1)
    xlabel!(L"x/a_x");ylabel!(L"y/a_y")
    title!(L"|\psi|^2")
    p2 = heatmap(x,y,angle.(ψ),aspectratio=1)
    xlabel!(L"x/a_x");ylabel!(L"y/a_y")
    title!(L"\textrm{phase}(\psi)")
    p = plot(p1,p2,size=(600,300))
    return p
end

# ==== set simulation parameters ====
sim = Par()
Lx = 20.0; Nx = 128
Ly = 20.0; Ny = 128
@pack! sim = Lx,Ly,Nx,Ny
initsim!(sim)
@unpack_Par sim
# ===================================

# declare the potential function
import FourierGPE.V
V(x,y,t)::Float64 = 0.5*(x^2 + y^2)

# useful TF state
ψ0(x,y,μ,g) = sqrt(μ/g)*sqrt(max(1.0-V(x,y,0.0)/μ,0.0)+im*0.0)
healing(x,y,μ,g) = 1/sqrt(g*abs2(ψ0(x,y,μ,g)))

#make convenient arrays
x,y,kx,ky,dx,dy,dkx,dky = makearrays(Lx,Nx,Ly,Ny)

#make initial state
ψi = ψ0.(x,y',μ,g)
ψi .+= (randn(Nx,Ny) |> complex)

# ====== Evolve in k space ==========
ϕi = kspace(ψi,sim)
sol = runsim(ϕi,sim)
# ===================================

# pull out the ground state:
ϕg = sol[end]
ψg = xspace(ϕg,sim)
showpsi(ψg,x,y)

# Add a vortex off-axis
using VortexDistributions

# initial state: imprint a vortex inside Thomas-Fermi radius
Rtf = sqrt(2*μ)
rv = 0.5*Rtf
xv,yv,cv = rv, 0.0, 1
initialvortex = [xv yv cv]
ξv = healing(xv,yv,μ,g)
ψv = copy(ψg)
makeallvortices!(ψv,initialvortex,x,y,ξv)
showpsi(ψv,x,y)

# In TF regim precession frequency is given analytically by:
# (see Fetter JLTP 2010)
ξ = 1/sqrt(μ)
Ωm = 3*log(Rtf/ξ/sqrt(2))/2/Rtf^2
Ωv = Ωm/(1-rv^2/Rtf^2)

#or a period of
Tv = 2*π/Ωv

# ==== set simulation parameters ====
γ = 0.0
tf = Tv
t = LinRange(ti,tf,Nt)
@pack! sim = tf,t,γ
initsim!(sim)
@unpack_Par sim
# ===================================

# ====== Evolve in k space ==========
ϕv = kspace(ψv,sim)
solv = runsim(ϕv,sim)
# ===================================

ϕf = solv[200]
ψf = xspace(ϕf,sim)
showpsi(ψf,x,y)

#trim last few frames to show one orbit
# analytical result is within 10%
anim = @animate for i=1:Nt-20
    showpsi(xspace(solv[i],sim),x,y)
end

gif(anim,"./examples/vortex.gif",fps=30)
