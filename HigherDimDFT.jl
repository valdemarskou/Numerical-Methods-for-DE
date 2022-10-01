using FFTW
using Plots
include("OtherDFTAlgorithms.jl")


# Forward 2 dimensional DFT. Requires both M and N to be even.

function Forward2DFFT(Array::{Float64,2})
    
end

A = Array{Float64,2}(rand(2^9,2^9))

fft(A)