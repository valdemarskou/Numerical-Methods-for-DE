


function DiscreteFourierTransform(fvalues::Array{Float64,1},s::Int)
    N = length(fvalues)
    FTransform = Array{ComplexF64,1}()
    for k=1:N 
        push!(FTransform,0)
        for j=1:N 
            FTransform[k] = FTransform[k] + fvalues[j] * exp(-2*s*pi*im*j*k/N)
        end
    end
    return FTransform

end


function InitializeFFT(N::Int,s::Int)
    w_first = exp(-2*s*pi*im/N)
    w = Array{ComplexF64,1}()
    for j=1:N 
        push!(w,(w_first)^j)
    end

    return w
end

# Anvend kun med værdier N=2^k hvor k er ulige.

function Radix2FFT(fvalues::Array{Float64,1},wvalues::Array{ComplexF64,1})
    # noPtsAtLevel,a,b,c,d,p,N,Ndiv2,m = repeat([Int],9)
    N = length(fvalues)
    Ndiv2 = Int(N/2)
    m = Int(log2(N))

    for l=1:Int((m+1)/2)
        noPtsAtLevel = Int(2^(l-1))
        a = 1
        b = Int(Ndiv2/noPtsAtLevel)
        return noPtsAtLevel
        p = 1
        for k=0:(b-1) 
            W = wvalues[p]
            for i=k:div(N, noPtsAtLevel):N
                z = W * (fvalues[a+i]-fvalues[b+i])
                fvalues[a+i] = fvalues[a+i] + fvalues[b+i]
                fvalues[b+i] = z
            end
            p = p + noPtsAtLevel
        end
    end
    #=
    for l=(div(m+3,2)):m 
        noPtsAtLevel = 2^(l-1)
        b = div(Ndiv2,noPtsAtLevel)
        c = noPtsAtLevel
        d = b + noPtsAtLevel
        p = 1

        for k=1:(b-1) 
            W = wvalues[p]
            for j=k:(N/noPtsAtLevel):(noPtsAtLevel - 1)
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
    =#
    return fvalues
    
end




test = [2.0,3,2,3,4,2,1,5]
# testOutput = DiscreteFourierTransform(V,1)

 wvalues = InitializeFFT(2^3,1)

Radix2FFT(test,wvalues)


