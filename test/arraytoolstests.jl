using RecursiveArrayTools
y = randn(3,3)
z = randn(4,4)
w = randn(5,5)+im*randn(5,5)
A = ArrayPartition(y,z)

# copyat_or_push!{T}(a::AbstractVector{T},i::Int,x)
# If i<length(x), it's simply a recursivecopy! to the ith element. Otherwise it will push! a deepcopy.
length(A)
length(y)+length()
copyat_or_push!(y,length(A),A)
