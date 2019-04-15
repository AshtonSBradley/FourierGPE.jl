function test()
  y = 3.
  x = 0.
  mysim = (:x,:y)
  myout = (:T,:S)
  # @eval $Tij = $dμi*$plantrans(ψtest,flags=flags)
for (out,sim) ∈ zip(myout,mysim)
    @eval $out = $sim
end
end

test()

y = 3.
x = 0.
mysim = (:x,:y)
myout = (:T,:S)
# @eval $Tij = $dμi*$plantrans(ψtest,flags=flags)
for (out,sim) ∈ zip(myout,mysim)
  @eval $out = $sim
end
