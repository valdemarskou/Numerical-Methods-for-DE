
# Algorithm for computing the values f(x_j) where x_j = 2*pi*j/N 
# for j=0..N-1, given the function f and the number of values N.

function AlmostEqual(a,b)
    if abs(a-b) < 0.0000000001
        return true
    else return false
    end 
end




function ComputeNodes(N::Int)
    xvalues = Array{Float64,1}(zeros(N))
    for j in 0:(N-1)
        xvalues[j+1] = 2*pi*j/N
    end
    return xvalues
end

function ComputeValues(N::Int, f::Function)
    fvalues = Array{Float64,1}(fill(0,N))
    for j=0:(N-1) 
        fvalues[j+1] = f(2*pi*j/N)
    end

    return fvalues
end

# Algorithm for directly evaluating the discrete fourier coefficients of f,
# given the values f(x_j) where x_j = 2*pi*j/N for j=0..N-1.

function DiscreteFourierCoefficients(fvalues::Array{Float64,1})
    N = length(fvalues)
    Ndiv2 = Int(N/2)
    fcoefficients = Array{ComplexF64,1}(fill(0,N+1))    
    for k in (-Ndiv2):Ndiv2
        s=0
        for j in 0:(N-1) 
            s=s+fvalues[j+1] * exp(-2*pi*im*j*k/N)
        end
        fcoefficients[k+Ndiv2+1] = s/N
    end

    return fcoefficients
end

# Algorithm for computing the fourier interpolant as a function
# of x, given the modes of f (computed through earlier algorithms)

function FourierInterpolantFromModes(x,fmodes::Array{ComplexF64,1})
    N = length(fmodes)-1
    # N will be even.
    Ndiv2 = Int(N/2)
    function f(x)
        return (fmodes[1]*exp(-im*x*Ndiv2) + fmodes[N+1]*exp(im*x*Ndiv2))/2
    end
    #=
    for k in (1-Ndiv2):(Ndiv2-1) 
        f = f + fmodes[k + 1 + Ndiv2] * exp(im*k*x)
    end
    =#
    g = Real(f)
    return g
end

# Algorithm for computing the fourier interpolant as a function
# of x, given the nodes of f (computed through earlier algorithms)

function FourierInterpolantFromNodes(x,xvalues,fvalues)
    N = length(xvalues)

    for j in 0:(N-1)
        if AlmostEqual(x,xvalues[j+1]) == true
            return fvalues[j+1]
        end      
    end
    s = 0
    for j in 0:(N-1)
        t = (x-xvalues[j+1])/2
        s = s+fvalues[j+1]*sin(N*t)*cot(t)/N
    end
    return s
end


