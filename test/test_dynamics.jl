# dynamics tests
tf = 10.0
t = LinRange(0.,tf,Nt)
sim2 = Sim(sim;γ = 0.0,tf=tf,t=t)

#make initial state
ϕi = sol[end]
@pack! sim2 = ϕi

sol2,err = testsim(sim2)
@test err == false

psitest = xspace(sol2[end],sim2)
npeak = abs2.(psitest[256])
@test isapprox(g*npeak,μ,rtol=1e-2)
