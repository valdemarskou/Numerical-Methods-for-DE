
using DifferentialEquations

function DiscreteFourierTransform(fvalues::Array{ComplexF64,1},s::Int)
    N = length(fvalues)
    DFT = Array{ComplexF64,1}()
    for k=1:N 
        push!(DFT,0)
        for j=1:N 
            DFT[k] = DFT[k] + fvalues[j] * exp(-2*s*pi*im*(j-1)*(k-1)/N)
        end
    end
    return DFT

end


function InitializeFFT(N::Int,s::Int)
    w_first = exp(-2*s*pi*im/N)
    w = Array{ComplexF64,1}()
    for j=1:N 
        push!(w,(w_first)^(j-1))
    end

    return w
end

# Anvend kun med værdier N=2^k hvor k er ulige.

function Radix2FFT(fvalues::Array{ComplexF64,1},wvalues::Array{ComplexF64,1})
    # noPtsAtLevel,a,b,c,d,p,N,Ndiv2,m = repeat([Int],9)
    N = length(fvalues)
    Ndiv2 = Int(N/2)
    m = Int(log2(N))

    for l=1:Int((m+1)/2)
        noPtsAtLevel = Int(2^(l-1))
        a = 0
        b = Int(Ndiv2/noPtsAtLevel)
        # return noPtsAtLevel
        p = 1
        for k=1:b 
            W = wvalues[p]
            for i=k:Int(N/noPtsAtLevel):(N-1)
                z = W * (fvalues[a+i]-fvalues[b+i])
                fvalues[a+i] = fvalues[a+i] + fvalues[b+i]
                fvalues[b+i] = z
            end
            p = p + noPtsAtLevel
        end
    end
    
    for l=Int((m+3)/2):m 
        noPtsAtLevel = Int(2^(l-1))
        b = Int(Ndiv2/noPtsAtLevel)
        c = noPtsAtLevel
        d = b + noPtsAtLevel
        p = 1
        a = 0

        for k=1:b 
            W = wvalues[p]
            for j=k:Int(N/noPtsAtLevel):(noPtsAtLevel - 1)
                for i=j:(2*noPtsAtLevel):N 
                    z=W*(fvalues[a+i]-fvalues[b+i])
                    fvalues[a+i] = fvalues[a+i] + fvalues[b+i]
                    fvalues[b+i] = fvalues[c+i] + fvalues[d+i]
                    fvalues[d+i] = W*(fvalues[c+i]-fvalues[d+i])
                    fvalues[c+i] = z
                end 
            end
            p = p+noPtsAtLevel
        end
    end
    
    return fvalues
    
end




test = Array{ComplexF64,1}([2.0,3.0,2,3.0,4,2,1,5])



# testOutput = DiscreteFourierTransform(test,1)

wvalues = InitializeFFT(2^8,1)

# testOutput2 = Radix2FFT(test,wvalues)




# Burde være nul.....
#testOutput[7]


testInput = Array{ComplexF64,1}(sin.(1:2^9)+cos.(1:2^9))


testOutput = Radix2FFT(testInput,wvalues)


testOutput[2]


# using Plots

#scatter(Real(testOutput))