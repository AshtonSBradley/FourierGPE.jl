# evolution
function V(x,y,t)::Float64
    return 0.5*(x^2 + y^2)*0.01
end

# nonlinearity
function nlin(ϕ,sim,t)
    @unpack g,x,y = sim
    ψ = xspace(ϕ,sim)
    ψ .*= @. g*abs2(ψ) + V(x,y',t)
    return kspace(ψ,sim)
end

function nlin!(dϕ,ϕ,sim,t)
    @unpack g,x,y = sim
    dϕ .= ϕ
    xspace!(dϕ,sim)
    dϕ .*= @. g*abs2(dϕ) + V(x,y',t)
    kspace!(dϕ,sim)
    return nothing
end

function Lgp!(dϕ,ϕ,sim,t)
    @unpack μ,γ,k2 = sim
    nlin!(dϕ,ϕ,sim,t)
    dϕ .+= @. (0.5*k2 - μ)*ϕ
    dϕ .*= @. -im*(1.0 - im*γ)
    return nothing
end

function Lgp(ϕ,sim,t)
    @unpack μ,γ,k2 = sim
    chi = nlin(ϕ,sim,t)
    return @. -im*(1.0 -im *γ)*((0.5*k2 - μ)*ϕ + chi )
end

function internalnorm(u,t)
    return sum((abs2.(u) .> 1e-6*maximum(abs2.(u))).*abs2.(u))
end

function runsim(ϕ,sim)
    prob = ODEProblem(Lgp!,ϕ,(sim.ti,sim.tf),sim)
    @info "𝒅𝜳 ⭆ Evolving in kspace"
    @info "damping γ = $(sim.γ)"
    @time sol = solve(prob,alg=Tsit5(),saveat=sim.t,reltol=1e-7)
    @info "⭆ Finished."
return sol
end
