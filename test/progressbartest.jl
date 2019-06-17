using OrdinaryDiffEq, Plots, LaTeXStrings, ParameterizedFunctions, Interact
gr(size=(600,200),grid=false);

#system parameters
λ = .7   # 🐇 linear growth
γ = 0.4  # 🐇 nonlinear loss
κ = 0.3  # 🐺 nonlinear growth
β = 1.2  # 🐺 linear loss

# 🎳 cute definition
pprey = @ode_def_all LotkaVolterra begin
  dr = λ*r - γ*w*r
  dw = - β*w + κ*r*w
end λ γ β κ

# here we declare what we want this function to depend upon parametrically
# we also need to package up all the parameters for later
params = [λ γ β κ]'

# ⚙️ in detail: (doing it without ParameterizedFunctions)
# function pprey(du,u,p,t)
#     du[1] =  λ*u[1] - γ*u[2]*u[1]
#     du[2] = -β*u[2] + κ*u[1]*u[2]
# end


#If we start at a fixed point, it is really fixed?
r̄ = β/κ+.5
w̄ = λ/γ-.7
u0 = [r̄ ; w̄]
ti = 0.0
tf = 2000000
Nt = 100
tspan = (ti,tf)
t = LinRange(ti,tf,Nt)
alg = Tsit5()

#define the problem type
prob = ODEProblem(pprey,u0,tspan,params)

#solve the problem, saving at specified time points
sol = solve(prob,alg=alg,saveat=t,progress=true)
