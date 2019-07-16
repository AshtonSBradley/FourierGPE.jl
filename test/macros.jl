


macro fun1(exi)
    ex = :( g(x) = $exi )
    return ex
end

macro test(Vex)
    ex = :( V(x,y,z,t)=$Vex )
    @show ex.args
end

@test y^4
V(0.1,0.1,0,0)

macro pot(exi)
    ex = :(V(x,y,z,t)=$exi)
    return ex
end

@pot x^2
V.(.1,0,0,0)


using Parameters

x=gensym()

eval(:($(x) = @with_kw (x=2,)))

getfield(Main,x)()

x




@spec(2)

ex = :( S= x1^2 + x2^2)
eval(ex(.1,.2))
typeof(ex)


function ksyms(n)
    ksyms = Symbol[]
    for j=1:n
        push!(ksyms,Symbol("k",j))
    end
return ksyms
end

ks = ksyms(20)

ex = :( x^2 for j=1:3  )
eval(ex)
# macro @kspec(n)
#     ks = ksyms(n)
#     ex = :([k^2,+) for k in ksyms )


    ex_new = ex
quote
    still_expression = $(esc(ex_new))
end


# You cannot "get a character" from a symbol
# You can use string to change a symbol to a string
# You can parse a string to a number
# Check out the args field of an expression: it's an array that holds a bunch of goodies!

k = ((1,2),(3,4))
ax = [kx^2 for kx in k, kx in k]

map(x -> x + 1, [1, 2, 3])

[x + 1 for x in [1, 2, 3]]

kvals = ((1.,2.,3.),(4.,5.,6.))


kx = kvec(1.,100)


# =====================================
# from D P Sanders intermediate julia
# =====================================

function wilkinson(n, x)
    result = x - 1

    for i in 2:n
        result *= (x - i)
    end

    result
end

p(n) = x -> wilkinson(n, x)

p(10)(.1)

ex = quote
        (x - 1) * (x - 2)
    end

ex = :( (x-1) * (x-2) )

typeof(ex)

ex

eval(ex)

x = 3.5
(x-1) * (x-2)

formula = "3x^2 + 4x - 3sin(x)"

eval(formula)

formula2 = Meta.parse(formula)


ex = :(3x^2 + 4x - 3sin(x))

typeof(ex)

eval(ex)

dump(ex)

Meta.show_sexpr(ex)

ex = :( (x-1) * (x-2) )

dump(ex)

ex.args[2]

ex2 = :(a=b+c)

ex2.head

ex.args[2].args[2]

typeof(ex.args[2].args[2])

ex2 = copy(ex)

ex2.args[2] = :z

ex2

ex

ex2.args[3].args[2]

ex2.args[3].args[2]=:z

ex2

z=1.2;x=3.2

eval(ex)

eval(ex2)

# Walking a syntax tree

# iterators product solution
things = [ [1,2], [10, 20, 30] ]
sum.(Iterators.product(things...))

things2 = ( (1,2), (10, 20, 30) ,(2.3,4.5,6.,7.))

kvals = Iterators.product(things...) |> collect

abs2.((1,2))

kind = Iterators.product(things...)
f(x) = sum(abs2.(x))
map(f,kind)
Iterators.sum(kind)

xvec(L,N)=LinRange(-L/2,L/2,N+1)[2:end] |> collect

function xvecs(L,N)
    X = []
    for (λ,ν) in zip(L,N)
        x = xvec(λ,ν)
        push!(X,x)
    end
    return X |> Tuple
end

function k2(L,N)
    M = length(L)
    @show M
    K = kvecs(L,N)
    if M==1
        K2 = [k^2 for k ∈ K[1] ]
    elseif M==2
        K2 = [kx^2 + ky^2 for kx in K[1], ky in K[2]]
    elseif M==3
        K2 = [kx^2 + ky^2 + kz^2 for kx in K[1], ky in K[2], kz in K[3]]
    elseif M==4
        K2 = [kx^2 + ky^2 + kz^2 + kw^2 for kx in K[1], ky in K[2], kz in K[3], kw in K[4]]
    end
    return K2 |> complex
end

function kvecs(L,N)
    K = []
    for (λ,ν) in zip(L,N)
        k = kvec(λ,ν)
        push!(K,k)
    end
    return K |> Tuple
end

L = (11.,)
N = (20,)
K = kvecs(L,N)

kind = Iterators.product(K...)
k3 = map(x-> sum(abs2.(x)),kind)

k1 = k2(L,N) |> real

all(k1==k3)


# macro for potential definition
# let's say we want to define the potential
#V(x,y,z,t) = 0.5*(x^2 + y^2 + 4*z^2)

using Revise, FourierGPE


import FourierGPE.V

V(x,y,z,t) = 0.5*(x^2 + y^2 + 10*z^2)

ex1 = :(V(x,y,z,t) = 0.5*(x^2 + y^2 + 10*z^2))

ex = :(
@with_kw mutable struct Params <: UserParams
    Vdef::Expr = :($ex1)
end
)

eval(ex)

par = Params()
eval(par.Vdef)
par.Vdef.args
V(.1,.1,.1,.1)


# what about build user parameters using macro?

macro params(Pname::Symbol,q)
    ex = :(q = $q)
    ex = :(
    @with_kw mutable struct $Pname <: UserParams @deftype Float64
        $ex
    end
    )
    return eval(ex)
end

a = .1
@params TestP5 a
par = TestP5()

@params mysim a
par = mysim()
@defparams TestP2 a
par = TestP()


# ===== easy start macro =====

@with_kw mutable struct Params <: UserParams @deftype Float64
    # user parameters:
    V::Expr = :( V(x,y,z,t) = 0.5*(x^2 + y^2 + 4*z^2) )
end
par = Params()

import FourierGPE.V
eval(par.V)

# ==== set simulation parameters ====
L=(15.,15.,15.)
N=(64,64,64)
γ = 0.05
tf = 1.5/γ
Nt = 200
t = LinRange(0.,tf,Nt)
# ========= Initialize simulation ======
sim = Sim(L,N,par)
@pack! sim = γ,tf,Nt,t
@unpack_Sim sim
