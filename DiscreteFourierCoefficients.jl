
# Algorithm for computing the values f(x_j) where x_j = 2*pi*j/N 
# for j=0..N-1, given the function f and the number of values N.

function ComputeValues(N::Int, f::Function)
    fvalues = Array{Float64,1}()
    for j=1:N 
        push!(fvalues,f(2*pi*j/N))
    end

    return fvalues
end

# Algorithm for directly evaluating the discrete fourier coefficients of f,
# given the values f(x_j) where x_j = 2*pi*j/N for j=0..N-1.

function DiscreteFourierCoefficients(fvalues::Array{Float64,1})
    fcoefficients = Array{ComplexF64,1}()    
    for k=(-length(fvalues)/2):(length(fvalues)/2)
        s=0
        for j=1:(length(fvalues)) 
            s=s+fvalues[j] * exp(-2*pi*im*j*k/(length(fvalues)))
        end
        push!(fcoefficients,s/(length(fvalues)))
    end

    return fcoefficients
end

# Algorithm for computing the fourier interpolant as a function
# of x, given the modes of f (computed through earlier algorithms)

function FourierInterpolantFromNodes(x,fmodes::Array{ComplexF64,1})
    N = length(fmodes)
    function f(x)
        return (fmodes[1]*exp(-im*N*x/2) + fmodes[N]*exp(im*N*x/2))/2
    end
    
    for k=2:(N-1) 
        f = f + fmodes[k]*exp(im*k*x)
    end
    g = Real(f)
    return g
end


# Testing env. with some function, som dog ikke er 2pi-periodisk

function f(x)
    return x*x+ 2*x
end


fvalues = ComputeValues(4,f)
DiscreteFourierCoefficients(fvalues)

