using Pkg, Revise

using FourierGPE

# ==== Units: ========================
# Oscillator

# ==== define user parameters =======
@with_kw mutable struct Params <: UserParams @deftype Float64
    # user parameters:
    V::Expr = :( V(x,y,z,t) = 0.5*(x^2 + y^2 + z^2) )
end
par = Params()

import FourierGPE.V
eval(par.V)

# ==== set simulation parameters ====
L = (15.,15.,15.)
N = (64,64,64)
γ = 0.5
tf = 2/γ
Nt = 200
t = LinRange(0.,tf,Nt)

# ========= Initialize simulation ======
# flags = FFTW.PATIENT
sim = Sim(L,N,par)
# sim = Sim{length(L)}(L=L,N=N,flags=flags,params=par)
@pack! sim = γ,tf,Nt,t
@unpack_Sim sim
sim

# ========= useful state functions ========
ψ0(x,y,z,μ,g) = sqrt(μ/g)*sqrt(max(1.0-V(x,y,z,0.0)/μ,0.0)+im*0.0)
healing(x,y,z,μ,g) = 1/sqrt(g*abs2(ψ0(x,y,z,μ,g)))

x,y,z = X
#make initial state
ψi = ψ0.(x,y',reshape(z,1,1,length(z)),μ,g)
ψi = randn(N)+im*randn(N)
ϕi = kspace(ψi,sim)
@pack! sim = ϕi,γ,tf,t

# ====== Evolve in k space ==========
sol = runsim(sim)
# ===================================


# ====== show slice using Plots ===
gr(titlefontsize=12,size=(500,300),colorbar=false)

# pull out the ground state
ϕg = sol[end]
ψg = xspace(ϕg,sim)
zmid = findfirst(z .≈ 0)
showpsi(x,y,ψg[:,:,zmid])

plot(x,y,abs2.(ψg[:,:,zmid]),linetype=:surface)

# === Dynamics of Kohn mode =======================
# kick it

ki = .7
ϕk = copy(ϕg)
ψk = xspace(ϕk,sim)
ψk = ψk.*exp.(im*ki*x)
ϕk = kspace(ψk,sim)
ϕi = copy(ϕk)
γ = 0.0
tf = 2*pi
t = LinRange(0.,tf,Nt)

simk = Sim(L,N,par)
@pack! simk = ϕi,γ,tf,Nt,t
@unpack_Sim simk

@time solk = runsim(simk)

# animate a slice
anim = @animate for i=1:Nt
    ϕ = solk[i]
    zmid = findfirst(z .≈ 0)
    ψ = xspace(ϕ,simk)[:,:,zmid]
    showpsi(x,y,ψ)
end

gif(anim,"./examples/3dKohn.gif",fps=30)
# gif(anim,"./examples/3dKohn.mov",fps=30)

# adimate slice as a surface
anim = @animate for i=1:Nt
    ϕ = solk[i]
    zmid = findfirst(z .≈ 0)
    ψ = xspace(ϕ,simk)[:,:,zmid]
    plot(x,y,abs2.(ψ),linetype=:surface)
end

gif(anim,"./examples/3dKohnSurface.mov",fps=30)

simk

# ========== animate isosurface in Makie ========
using Makie, AbstractPlotting

function dense(phi)
    ψm = xspace(phi,sim)
    density = abs2.(ψm)
    pmax = maximum(density)
    return density/pmax
end

function densityfilm(sol,Nt)
    saveto="examples/3dKohn.gif"
    scene = Scene()
    tindex = Node(1)
    scene = volume(lift(i -> dense(sol[i]), tindex),
    algorithm = :iso,
    show_axis=false,
    isovalue=3f0(.15))

    R = 70
    eyepos = Vec3f0(R,R,R)
    lookat = Vec3f0(18,18,0)

    record(scene, saveto, 1:Nt) do i
        update_cam!(scene, eyepos, lookat)
        rotate_cam!(scene, 0., -0.25, 0.)
        tindex[] = i
    end
    p = scene[end];
    return
end

p = densityfilm(sol,Nt)

# p.attributes.attribute
