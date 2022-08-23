

include("OtherDFTAlgorithms.jl")


# Evaluating the k'th Legendre Polynomial at a real value x. Unfortunately, I cannot
# figure out how to return the abstract function, at least for the moment.

function LegendrePolynomial(k::Int, x)
    if k==0
        return 1
    end

    if k==1
        return x
    end

    Lkmin2 = 1
    Lkmin1 = x
    Lkmin0 = Float64

    for j in 2:k 
      Lkmin0 = x*Lkmin1*(2*j-1)/j - Lkmin2 * (j-1)/j   
        Lkmin2 = Lkmin1
        Lkmin1 = Lkmin0
    end

    return Lkmin0
end

# Evaluating the k'th Chebyshev Polynomial at a real value x. The cutoff value for the  Unfortunately, I cannot
# figure out how to return the abstract function, at least for the moment.

function ChebyshevPolynomial(k::Int, x)
    if k == 0
        return 1
    end

    if k == 1
        return x
    end

    Tkmin0 = Float64

    if k < 70
        Tkmin2 = 1
        Tkmin1 = x
        for j in 2:k
            Tkmin0 = 2*x*Tkmin1 - Tkmin2
            Tkmin2 = Tkmin1
            Tkmin1 = Tkmin0 
        end
    
    else
        Tkmin0 = cos(k*acos(x))
    end

    return Tkmin0
end


# Evaluating the k'th Legendre Polynomial and its derivative at a real value x. 

function LegendrePolynomialAndDerivative(N::Int, x)

    LN = Float64
    LNdot = Float64

    if N == 0
        LN = 1
        LNdot = 0

    elseif N == 1
        LN = x
        LNdot = 1

    else
        LNmin2 =  1
        LNmin1 = x
        LNdotmin2 = 0
        LNdotmin1 = 1

        for k in 2:N
            LN = x * LNmin1 * (2k-1)/(k) - LNmin2 * (k-1)/(k)
            LNdot = LNdotmin2 + (2k-1) * LNmin1
            LNmin2 = LNmin1
            LNmin1 = LN
            LNdotmin2 = LNdotmin1
            LNdotmin1 = LNdot 
        end
    
    end

    return LN, LNdot

end


# Evaluating Legendre Gauss Nodes and Weights. We set n_it = 4, as well
# as TOL = 4*epsilon. Hopefully it works, not sure how to check it in practice.

function LegendreGaussNodesAndWeights(N::Int,epsilon::Float64)
    xvalues = Array{Float64,1}(zeros(N+1))
    wvalues = Array{Float64,1}(zeros(N+1))

    if N == 0
        xvalues[1] = 0
        wvalues[1] = 2
    
    elseif N == 1
        xvalues[1] = -sqrt(1/3)
        wvalues[1] = 1
        xvalues[2] = -xvalues[1]
        wvalues[2] = wvalues[1]

    else
        for j in 0:Int((floor((N+1)/2)-1)) 
            xvalues[j+1] = -cos((2j+1)/(2*N+2) * pi)

            for k in 1:4
                LNplus1, LdotNplus1 = LegendrePolynomialAndDerivative(N+1, xvalues[j+1])
                delta = - LNplus1/LdotNplus1
                xvalues[j+1] = xvalues[j+1] + delta
                if abs(delta) < 4*epsilon * abs(xvalues[j+1])
                    break
                end
            end
            LNplus1, LdotNplus1 = LegendrePolynomialAndDerivative(N+1, xvalues[j+1])
            xvalues[N+1-j] = -xvalues[j+1]
            wvalues[j+1] = 2/((1-xvalues[j+1]^2)*LdotNplus1^2)
            wvalues[N+1-j] = -wvalues[j+1]
        end
    end

    if mod(N, 2) == 0
        LNplus1, LdotNplus1 = LegendrePolynomialAndDerivative(N+1, 0.0)
        xvalues[Int(N/2)+1] = 0
        wvalues[Int(N/2)+1] = 2/(LdotNplus1^2)
    end

    return xvalues, wvalues
end


# Combined algorithm to compute L_N(x) and q(x) = L_n+1(x) - L_n-1(x) and q'(x).

function qAndLEvaluation(N::Int,x)
    LNmin2 = 1
    LNmin1 = x
    LN = Float64
    LdotNmin2 = 0
    LdotNmin1 = 1
    LdotN = Float64
    
    for k in 2:N 
        LN = x * LNmin1 * (2k-1)/(k) - LNmin2 * (k-1)/(k)
        LdotN = LdotNmin2 + (2k-1) * LNmin1
        LNmin2 = LNmin1
        LNmin1 = LN
        LdotNmin2 = LdotNmin1
        LdotNmin1 = LdotN
    end
    k = N+1
    LNplus1 = (2k-1)/(k) * x * LN - (k-1)/(k) * LNmin1
    LdotNplus1 = LdotNmin1 + (2k-1) * LN
    q = LNplus1 - LNmin1
    qdot = LdotNplus1 -LdotNmin1
    
    return q, qdot, LN
end


# Evaluating the Legendre Gauss-Lobatto Nodes and Weights. We again set n_it = 4, as
# well as TOL = 4*epsilon.

function LegendreGaussLobattoNodesAndWeights(N::Int,epsilon::Float64)
    xvalues = Array{Float64,1}(zeros(N+1))
    wvalues = Array{Float64,1}(zeros(N+1))
    
    if N == 1
        xvalues[1] = -1
        wvalues[1] = 1 
        xvalues[2] = 1
        wvalues[2] = wvalues[1]

    else  
        xvalues[1] = -1
        wvalues[1] = 2/(N*(N+1))
        xvalues[N+1] = 1
        wvalues[N+1] = wvalues[1]

        for j in 1:Int((floor((N+1)/2)-1)) 
            xvalues[j+1] = -cos((j+1/4)*pi/N - 3/((j+1/4)*8*N*pi))

            for k in 0:4 
                q,qdot,LN = qAndLEvaluation(N,xvalues[j+1])
                delta = -q/qdot
                xvalues[j+1] = xvalues[j+1] + delta
                if abs(delta) < 4*epsilon*abs(xvalues[j+1])
                    break
                end
            end
            q,qdot,LN = qAndLEvaluation(N,xvalues[j+1])
            xvalues[N+1-j] = -xvalues[j+1]
            wvalues[j+1] = 2/(N*(N+1)*LN^2)
            wvalues[N+1-j] = wvalues[j+1]
        end
    end
    if mod(N,2) == 0
        q,qdot,LN = qAndLEvaluation(N,0.0)
        xvalues[Int(N/2)+1] = 0
        wvalues[Int(N/2)+1] = 2/(N*(N+1)*LN^2)
    end

    return xvalues, wvalues
end


# Evaluating the Chebyshev Gauss Nodes and Weights.
#

function ChebyshevGaussNodesAndWeights(N::Int)
    xvalues = Array{Float64,1}(zeros(N+1))
    wvalues = Array{Float64,1}(zeros(N+1))

    for j in 0:N
        xvalues[j+1] = -cos(pi*(2j+1)/(2*N+2))
        wvalues[j+1] = pi/(N+1) 
    end
    
    return xvalues, wvalues
end


# Evaluating the Chebyshev Gauss Lobatto Nodes and Weights.
# Algorithm 27 in the book.

function ChebyshevGaussLobattoNodesAndWeights(N::Int)
    xvalues = Array{Float64,1}(zeros(N+1))
    wvalues = Array{Float64,1}(zeros(N+1))

    for j in 0:N
        xvalues[j+1] = -cos(j*pi/N)
        wvalues[j+1] = pi/N
    end

    wvalues[1] = wvalues[1]/2
    wvalues[N+1] = wvalues[N+1]/2
    
    return xvalues, wvalues
end


# Initialization algorithm for the fast cosine transform. Computes the C and S vectors,
# in that order.

function InitializeFCT(N::Int)
    C = Array{Float64,1}(zeros(N+1))
    S = Array{Float64,1}(zeros(N+1))

    for j in 0:N 
        C[j] = cos(pi*j/N)
        S[j] = sin(pi*j/N)
    end
    
    return C, S
end


# Fast cosine transform, using the fft method from the FFTW.jl package.
# Algorithm 28. Virker igen kun når N er lige.

function FastCosineTransform(fvalues::Array{ComplexF64,1},C,S,s::Int)
    N = length(fvalues)-1   
    evalues = Array{ComplexF64,1}(zeros(N+1))

    for j in 0:(N-1) 
        evalues[j+1] = (fvalues[j+1] - fvalues[N-j+1])/2 + S[j+1] * (fvalues[j+1] - fvalues[N-j+1])
    end

    abarvalues, bbarvalues = ForwardRealFFT(evalues)
    avalues = Array{Float64,1}(zeros(N+1))

    for k in 0:Int(N/2) 
        avalues[2k+1] = abarvalues[k+1]
    end

    avalues[2] = fvalues[1] - fvalues[N+1]
    
    for j in 1:(N-1) 
        avalues[2] = avalues[2] + 2*C[j+1]*fvalues[j+1]
    end

    avalues[2] = avalues[2]/N

    for k in 1:(Int(N/2)-1) 
       avalues[2k+2] = bbarvalues[k+1] - avalues[2k] 
    end

    if s == -1
        for k in 0:N 
            avalues[k+1] = N*avalues[k+1]/2
        end
    end

    return avalues
end


# Fast Chebyshev transform.
# Algorithm 29. Virker igen kun når N er lige.

function FastChebyshevTransform(fvalues::Array{ComplexF64},C,S,s::Int)
    N = length(fvalues)-1
    gvalues = fvalues

    if s == -1
        gvalues[1] = 2*gvalues[1]
        gvalues[N+1] = 2*gvalues[N+1]
    end
    
    avalues = FastCosineTransform(gvalues,C,S,s)

    if s == 1
        avalues[1] = avalues[1]/2
        avalues[N+1] = avalues[N+1]/2
    end

    return avalues
end


# Evaluating the BarycentricWeights given the nodes (xvalues)
# Algorithm 30.

function BarycentricWeights(xvalues::Array{ComplexF64,1})
    N = length(xvalues)-1

    wvalues = Array{ComplexF64,1}(fill(1,N+1))
    
    for j in 1:N
        for k in 0:(j-1) 
            wvalues[k+1] = wvalues[k+1]*(xvalues[k+1] - xvalues[j+1])
            wvalues[j+1] = wvalues[j+1]*(xvalues[j+1] - xvalues[k+1])
        end 
    end

    for j in 0:N
        wvalues[j+1] = 1/(wvalues[j+1]) 
    end

    return wvalues
end

isapprox(1,1.0000001)
# Lagrange interpolant from barycentric form.
# Algorithm 31.

function LagrangeInterpolation(x,xvalues::Array{ComplexF64,1},fvalues::Array{ComplexF64,1},wvalues::Array{ComplexF64,1})
    N = length(xvalues)-1
    numerator = 0
    denominator = 0

    for j in 0:N
        if isapprox(x,xvalues[j+1]) == true
            return fvalues[j+1]
        end
        t = wvalues[j+1]/(x-xvalues[j+1])
        numerator = numerator + t*fvalues[j+1]
        denominator = denominator + t
    end

    return numerator/denominator
end


# Matrix for interpolation between two sets of points. xvalues and wvalues are of 
# length N+1. xivalues are of length M+1. Algorithm 32.


function PolymialInterpolationMatrix(xvalues::Array{ComplexF64,1},wvalues::Array{ComplexF64,1},ξ::Array{ComplexF64,1})
    @assert(length(xvalues) == length(wvalues))
    N = length(xvalues) - 1
    M = length(ξ) - 1
    T = Array{ComplexF64,2}(fill(0,(M+1,N+1)))
    for k in 0:M 
        rowHasMatch = false
        for j in 0:N
            if isapprox(ξ[k+1],xvalues[j+1]) == true
                rowHasMatch = true
                T[k+1,j+1] = 1
            end
        end
        if rowHasMatch == false
            s = 0
            for j in 0:N
                t = wvalues[j+1]/(ξ[k+1] - xvalues[j+1])
                T[k+1,j+1] = t
                s = s + t
            end

            for j in 0:N 
                T[k+1,j+1] = T[k+1,j+1]/s
            end
        end
    end

    return T
end


# Interpolation Between Two Sets of Points by Matrix Multiplication. Not sure if there
# is a less manual way of doing this. Algorithm 33.

function InterpolateToNewPoints(T::Array{ComplexF64,2},f::Array{ComplexF64,1})
    M,N = size(T)
    @assert(N == length(f)) 
    finterp = Array{ComplexF64,1}(fill(0,M+1))

    for i in 0:M
        t = 0
        for j in 0:N 
            t = t + T[i+1,j+1]*f[j+1]
        end
        finterp[i+1] = t
    end

    return finterp
end


#
# Algorithm 34.

function LagrangeInterpolatingPolynomials(x,xvalues::Array{ComplexF64,1},wvalues::Array{ComplexF64,1})
    # The vectors xvalues and wvalues are of the same length N. Can put assert statement.
    N = length(xvalues) - 1
    ℓ = Array{ComplexF64,1}(fill(0,N+1))
    
    xMatchesNode = false
    for j in 0:N
        if isapprox(x,xvalues[j+1]) == true
            ℓ[j+1] = 1
            xMatchesNode = true
        end
    end

    if xMatchesNode == true
        return ℓ
    end

    s = 0
    for j in 0:N
        t = wvalues[j+1]/(x-xvalues[j+1])
        ℓ[j+1] = t
        s = s + t
    end
    for j in 0:N
        ℓ[j+1] = ℓ[j+1]/s
    end

    return ℓ
end







# # # TESTING ENVIRONMENT # # #

test1 = Array{ComplexF64,1}(randn(2^5+1))
test2 = Array{ComplexF64,1}(randn(2^6+1))
weighttest1 = BarycentricWeights(test1)
T = PolymialInterpolationMatrix(test1,weighttest1,test2)
size(T)

LagrangeInterpolatingPolynomials(3,test1,weighttest1)