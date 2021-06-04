#build global const transform lib as in:
# https://github.com/JuliaLang/julia/blob/8535b8c4d7f4179f7b13a04f16074a2446dda3ac/stdlib/Random/src/RNGs.jl#L364-L366
# struct _GLOBAL_RNG <: AbstractRNG
#     global const GLOBAL_RNG = _GLOBAL_RNG.instance
# end



abstract type AbstractLib end

mutable struct TransformLib <: AbstractLib 
    t1::Matrix{Float64}
    t2::Matrix{Float64}
    
    TransformLib(t1,t2) = new(t1,t2)
    TransformLib() = TransformLib(rand([0 1],2,2),rand([0 1],2,2))
    global const GLOBAL_TLIB = TransformLib()
end


##
