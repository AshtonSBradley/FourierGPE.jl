# initialize transforms and fields
function initsim!(sim;flags=FFTW.MEASURE)
    @unpack L,N = sim
    X,K,dX,dK,DX,DK,T = makearraystransforms(L,N)
    espec = 0.5*k2(L,N)
    @pack! sim = T,X,K,espec
    return nothing
end

# evolution: default potential
V(x,t)::Float64 = 0.5*x^2
V(x,y,t)::Float64 = 0.5*(x^2 + y^2)
V(x,y,z,t)::Float64 = 0.5*(x^2 + y^2 + z^2)

# nonlinearity
function nlin(ϕ,sim,t)
    @unpack g,x,y = sim
    ψ = xspace(ϕ,sim)
    ψ .*= @. g*abs2(ψ) + V(x,y',t)
    return kspace(ψ,sim)
end

function nlin!(dϕ,ϕ,sim::Sim{1},t)
    @unpack g,X = sim; x = X[1]
    dϕ .= ϕ
    xspace!(dϕ,sim)
    @. dϕ *= g*abs2(dϕ) + V(x,t)
    kspace!(dϕ,sim)
    return nothing
end

function nlin!(dϕ,ϕ,sim::Sim{2},t)
    @unpack g,X = sim; x,y = X
    dϕ .= ϕ
    xspace!(dϕ,sim)
    @. dϕ *= g*abs2(dϕ) + V(x,y',t)
    kspace!(dϕ,sim)
    return nothing
end

function nlin!(dϕ,ϕ,sim::Sim{3},t)
    @unpack g,X = sim; x,y,z = X
    y = y'; z = reshape(z,(1,1,length(z)))
    dϕ .= ϕ
    xspace!(dϕ,sim)
    @. dϕ *= g*abs2(dϕ) + V(x,y,z,t)
    kspace!(dϕ,sim)
    return nothing
end

function Lgp!(dϕ,ϕ,sim,t)
    @unpack γ,μ,espec = sim
    nlin!(dϕ,ϕ,sim,t)
    @. dϕ = -im*(1.0 - im*γ)*(dϕ + (espec - μ)*ϕ)
    return nothing
end

function Lgp(ϕ,sim,t)
    @unpack μ,γ,espec = sim
    chi = nlin(ϕ,sim,t)
    return @. -im*(1.0 - im *γ)*((espec - μ)*ϕ + chi )
end

function internalnorm(u,t)
    return sum((abs2.(u) .> 1e-6*maximum(abs2.(u))).*abs2.(u))
end

function runsim(sim,ϕ=sim.ϕi;info=true)
    prob = ODEProblem(Lgp!,ϕ,(sim.ti,sim.tf),sim)
    info && @info "𝒅𝜳 ⭆ Evolving in kspace"
    info && @info "damping γ = $(sim.γ)"
    sol = solve(prob,alg=sim.alg,saveat=sim.t,reltol=sim.reltol)
    info && @info "⭆ Finished."
return sol
end
