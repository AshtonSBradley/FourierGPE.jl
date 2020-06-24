## Tests for struct evolution
using Pkg;Pkg.activate(".")

## Fix types
"""
A = crandn_array(M)

Make placeholder `2x2x...` complex `randn()` array of `M` dimensions."""
function crandn_array(M)
    a = Int.(ones(M)).+1
    return randn(a...) |> complex
end

using FFTW

a = crandn_array(2)
p1 = plan_fft(a)
typeof(p1)