
using DifferentialEquations
using FFTW

function DiscreteFourierTransform(fvalues::Array{ComplexF64,1},s::Int)
    N = length(fvalues)
    DFT = Array{ComplexF64,1}(zeros(N))
    for k in 0:(N-1)
        for j in 0:(N-1) 
            DFT[k+1] = DFT[k+1] + fvalues[j+1] * exp(-2*s*pi*im*j*k/N)
        end
    end
    return DFT

end


function InitializeFFT(N::Int,s::Int)
    # w_first = exp(-2*s*pi*im/N)
    w = Array{ComplexF64,1}(zeros(N))
    for j=1:N 
        #push!(w,exp(-2*s*pi*j*im/N))
        w[j] = exp(-2*s*pi*j*im/N)
    end

    return w
end

# Can only be applied to lengths N=2^k where k is odd. This seems to be correct 
# in terms of what is written in the textbook, however the values disagree with 
# the naive implementation of the DFT, as well as with the methods from the FFTW 
# package. Unsure what to do about this issue.

function Radix2FFT(fvalues::Array{ComplexF64,1},wvalues::Array{ComplexF64,1})
    # noPtsAtLevel,a,b,c,d,p,N,Ndiv2,m = repeat([Int],9)
    N = length(fvalues)
    Ndiv2 = Int(N/2)
    m = Int(log2(N))

    for l=1:Int((m+1)/2)
        noPtsAtLevel = Int(2^(l-1))
        a = 0
        b = Int(Ndiv2/noPtsAtLevel)
        # return b
        p = 0
        for k=0:(b-1) 
            W = wvalues[p+1]
            # return W
            for i=k:Int(N/noPtsAtLevel):(N-1)
                z = W * (fvalues[a+i+1]-fvalues[b+i+1])
                # return (fvalues[a+i+1]-fvalues[b+i+1])
                fvalues[a+i+1] = fvalues[a+i+1] + fvalues[b+i+1]
                fvalues[b+i+1] = z 
                return fvalues[1]
            end
            p = p + noPtsAtLevel
        end
    end
    
    for l=Int((m+3)/2):m
        noPtsAtLevel = Int(2^(l-1))
        a = 0
        b = Int(Ndiv2/noPtsAtLevel)
        c = noPtsAtLevel
        d = b + noPtsAtLevel
        p = 0

        for k=0:(b-1) 
            W = wvalues[p+1]
            for j=k:Int(N/noPtsAtLevel):(noPtsAtLevel - 1)
                for i=j:(2*noPtsAtLevel):(N-1) 
                    z = W*(fvalues[a+i+1]-fvalues[b+i+1])
                    fvalues[a+i+1] = fvalues[a+i+1] + fvalues[b+i+1]
                    fvalues[b+i+1] = fvalues[c+i+1] + fvalues[d+i+1]
                    fvalues[d+i+1] = W*(fvalues[c+i+1]-fvalues[d+i+1])
                    fvalues[c+i+1] = z
                end 
            end
            p = p+noPtsAtLevel
        end
    end
    
    return fvalues
    
end


#=
# Example. Delete if necessary.
testInput = Array{ComplexF64,1}(sin.(1:2^9)+cos.(1:2^9)+(1:2^9)/(2^9))
DFTOutput = DiscreteFourierTransform(testInput,1)

wvalues = InitializeFFT(2^9,1)
FFTOutput = Radix2FFT(testInput,wvalues)

# FFTWOutput = fft(testInput)

using Plots
scatter(DFTOutput,xlims=(-100,100),ylims=(-300,300))
# scatter(FFTOutput,xlims=(-100,100),ylims=(-300,300))
# scatter(FFTWOutput)

=#

