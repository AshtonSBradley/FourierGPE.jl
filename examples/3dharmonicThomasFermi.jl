using Plots
gr(titlefontsize=12,size=(500,300),colorbar=false)


using FourierGPE

# Units: oscillator
L = (20.,20.,20.)
N = (64,64,64)
sim = Sim(L,N)
@unpack_Sim sim

import FourierGPE.V
V(x,y,z,t) = 0.5*(x^2 + y^2 + 4*z^2)

## set simulation parameters ====
γ = 0.5
μ = 15.0
tf = 0.3 # 1.5/γ
Nt = 150
t = LinRange(0.,tf,Nt)

## state functions
ψ0(x,y,z,μ,g) = sqrt(μ/g)*sqrt(max(1.0-V(x,y,z,0.0)/μ,0.0)+im*0.0)
healing(x,y,z,μ,g) = 1/sqrt(g*abs2(ψ0(x,y,z,μ,g)))

## initial state
x,y,z = X
y = y'; z = reshape(z,1,1,length(z))
ψi = ψ0.(x,y,z,μ,g)
ϕi = kspace(ψi,sim)

@pack_Sim! sim

plot(x,abs2.(ψ0.(x,0.,0.,μ,g)))
plot!(x,abs2.(ψi[:,32,32]))

## evolve in k space
sol = runsim(sim)

## ground state
ϕg = sol[end]
ψg = xspace(ϕg,sim)
showpsi(x,y',ψg[:,:,32])

R(w) = sqrt(2*μ/w^2)
R(1)
R(2)

plot(x,abs2.(ψ0.(x,0.,0.,μ,g)))
plot!(x,abs2.(ψg[:,32,32]))
μ/g

abs2.(ψg[32,32,32])

# animate a slice
anim = @animate for i=1:Nt
    ϕ = sol[i]
    zmid = findfirst(z .≈ 0)
    ψ = xspace(ϕ,sim)
    showpsi(x,y,ψ[:,:,zmid])
end
gif(anim,"./examples/3dquench.gif",fps=30)
#TODO debug
## animate isosurface in Makie
using Makie, AbstractPlotting

function dense(phi)
    ψm = xspace(phi,sim)
    density = abs2.(ψm)
    pmax = maximum(density)
    return density/pmax
end

function densityfilm(sol,Nt;file="3dquench.gif")
    saveto="examples/"*file
    scene = Scene()
    tindex = Node(1)
    scene = volume(lift(i -> dense(sol[i]), tindex),
    algorithm = :iso,
    color = (c3,0.25),
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

p = densityfilm(sol,Nt,file="3dquench1.gif")

# p.attributes.attribute
