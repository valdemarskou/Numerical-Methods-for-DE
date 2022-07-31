



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
    noPtsAtLevel,a,b,c,d,p,N,Ndiv2,m=Int
    N = length(fvalues)
    Ndiv2 = N/2
    m = log2(N)

    for l=1:((m+1)/2) 
        noPtsAtLevel = 2^(l-1)
        a = 0
        b = Ndiv2/noPtsAtLevel
        p = 1
        for k=0:(b-1) 
            W = wvalues[p]
            for  in 
                
            end
        end
    end

    
end


# InitializeFFT(2^8,1)

V=[3,4,5]

V[2] = 2
