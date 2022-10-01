

include("OtherDFTAlgorithms.jl")
using SpecialFunctions

# Evaluate the Legendre Coefficients of the Derivative of a Polynomial. Both 
# for real and complex values. Algorithm 4.

function LegendreDerivativeCoefficients(f̂::Array{ComplexF64,1})
    
    N = length(f̂) - 1
    f̅¹ = Array{ComplexF64,1}(fill(0,N+1))

    f̅¹[N+1] = 0
    f̅¹[N] = (2*N - 1)*f̂[N+1]

    for k in reverse(0:(N-2))
        f̅¹[k+1] = (2*k + 1)*(f̂[k+2] + f̂[k+3]/(2*k+5))
    end

    return f̅¹
end

function LegendreDerivativeCoefficients(f̂::Array{Float64,1})
    
    N = length(f̂) - 1
    f̅¹ = Array{ComplexF64,1}(fill(0,N+1))

    f̅¹[N+1] = 0
    f̅¹[N] = (2*N - 1)*f̂[N+1]

    for k in reverse(0:(N-2))
        f̅¹[k+1] = (2*k + 1)*(f̂[k+2] + f̂[k+3]/(2*k+5))
    end

    return f̅¹
end



#Evaluate the Chebyshev Coefficients of the Derivative of a Polynomial.
# Algorithm 5.

function ChebyshevDerivativeCoefficients(f̂::Array{ComplexF64,1})

    N = length(f̂) - 1
    f̅¹ = Array{ComplexF64,1}(fill(0,N+1))

    f̅¹[N+1] = 0
    f̅¹[N] = 2*N*f̂[N+1]

    for k in reverse(0:(N-2))
        f̅¹[k+1] = 2*(k+1)*f̂[k+2] + f̅¹[k+3]       
    end

    f̅¹[1] = f̂[2] + f̅¹[3]/2

    return f̅¹
end

function ChebyshevDerivativeCoefficients(f̂::Array{Float64,1})

    N = length(f̂) - 1
    f̅¹ = Array{ComplexF64,1}(fill(0,N+1))

    f̅¹[N+1] = 0
    f̅¹[N] = 2*N*f̂[N+1]

    for k in reverse(0:(N-2))
        f̅¹[k+1] = 2*(k+1)*f̂[k+2] + f̅¹[k+3]       
    end

    f̅¹[1] = f̂[2] + f̅¹[3]/2

    return f̅¹
end

# Unnormalized I think. Finally works. 
#

function JacobiPolynomial(a::Float64, b::Float64, k::Int, x)
    if(k == 0)
        return 1
    end

    if(k == 1)
        return 1/2 * (a - b + (a + b + 2)*x)
    end

    Pkmin2 = 1
    Pkmin1 = 1/2 * (a - b + (a + b + 2)*x)
    Pkmin0 = Float64

    for j in 2:k
        
        Pkmin0 = ((a*a-b*b)/((2*j+a+b)*(2*j+a+b-2))+x)*Pkmin1 - 2*(j-1+a)*(j-1+b)/((2*j+a+b-1)*(2*j+a+b-2))*Pkmin2
        Pkmin0 = Pkmin0 * (2*j+a+b)*(2*j+a+b-1)/(2*j*(j+a+b))
        Pkmin2 = Pkmin1
        Pkmin1 = Pkmin0
    end
    
    return Pkmin0
end

function OrthonormalJacobiPolynomial(a::Float64, b::Float64, k::Int, x)
    y = 2^(a+b+1) * (gamma(k+a+1)*gamma(k+b+1))/(gamma(k+1)*gamma(k+a+b+1)*(2*k+a+b+1))
    
    return JacobiPolynomial(a,b,k,x)/(sqrt(y))
end


# Derivative of the k'th Jacobi Polynomial. Based on the previous routine. Works.
#

function JacobiPolynomialDerivative(a::Float64, b::Float64, k::Int, x)
    if k==0
        return 0
    else
        return 1/2 * (a+b+k+1)*JacobiPolynomial(a+1,b+1,k-1,x)
    end
end

function OrthonormalJacobiPolynomialDerivative(a::Float64, b::Float64, k::Int, x)
    if k==0
        return 0
    else
        return sqrt(k*(k+a+b+1))*OrthonormalJacobiPolynomial(a+1,b+1,k-1,x)
    end

end


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
# 

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
# as TOL = 4*epsilon. Hopefully it works.

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
            xvalues[j+1] = -cos((2*j+1)/(2*N+2) * pi)

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


function DiscreteLegendreCoefficients(M::Int,xvalues::Array{Float64,1},wvalues::Array{Float64,1},f::Function)
    N = length(xvalues) - 1
    fvalues = f.(xvalues)
    fcoeffs = Array{Float64,1}(zeros(M))


    for k in 0:(M-1)
        numerator = 0
        #denominator = 0
        for j in 0:N
            numerator = numerator + fvalues[j+1]*LegendrePolynomial(k,xvalues[j+1])*wvalues[j+1]
            #denominator = denominator + wvalues[j+1]*LegendrePolynomial(k,xvalues[j+1])^2
        end
        fcoeffs[k+1] = numerator * (2*k+1)/2
    end
    
    return fcoeffs
end


# Initialization algorithm for the fast cosine transform. Computes the C and S vectors,
# in that order.

function InitializeFCT(N::Int)
    C = Array{Float64,1}(zeros(N+1))
    S = Array{Float64,1}(zeros(N+1))

    for j in 0:N 
        C[j+1] = cos(pi*j/N)
        S[j+1] = sin(pi*j/N)
    end
    
    return C, S
end


# Fast cosine transform, using the fft method from the FFTW.jl package.
# Algorithm 28. Virker igen kun når N er lige.

function FastCosineTransform(fvalues::Array{ComplexF64,1},C,S,s::Int)
    N = length(fvalues)-1   
    evalues = Array{Float64,1}(zeros(N))

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

function FastChebyshevTransform(fvalues::Array{ComplexF64,1},C,S,s::Int)
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

# Evaluating the Vandermonde matrix associated to (unnormalized) Jacobi polynomials 
# defined by a,b, and quadrature nodes xvalues. This might not be correct due to the
# unnormalization, but it should work for the Legendre polynomials.

function VandermondeMatrix(a::Float64,b::Float64,xvalues::Array{Float64,1})
    N = length(xvalues) - 1
    V = Array{Float64,2}(fill(0,N+1,N+1))

    for i in 0:N
        for j in 0:N
            #V[i+1,j+1] = JacobiPolynomial(a,b,j,xvalues[i+1])
            V[i+1,j+1] = OrthonormalJacobiPolynomial(a,b,j,xvalues[i+1])
        end
    end
    return V
end

function VandermondeDerivativeMatrix(a::Float64,b::Float64,xvalues::Array{Float64,1})
    N = length(xvalues) - 1
    V = Array{Float64,2}(fill(0,N+1,N+1))

    for i in 0:N
        for j in 0:N
            #V[i+1,j+1] = JacobiPolynomialDerivative(a,b,j,xvalues[i+1])
            V[i+1,j+1] = OrthonormalJacobiPolynomialDerivative(a,b,j,xvalues[i+1])
        end
    end
    return V
end


# Lagrange interpolant from the Vandermonde Matrix formula. Evaluates all N 
# interpolating polynomials. Works for Legendre polynomials.

function LagrangeInterpolant(a::Float64,b::Float64,xvalues::Array{Float64,1},x)
    N = length(xvalues) - 1
    V = inv(VandermondeMatrix(a,b,xvalues))
    h = Array{Float64,1}(zeros(N+1))

    for j in 0:N
        t = 0
        for n in 0:N
            #t = t+V[n+1,j+1]*JacobiPolynomial(a,b,n,x)
            t = t+V[n+1,j+1]*OrthonormalJacobiPolynomial(a,b,n,x)
        end
        h[j+1] = t
    end

    return h
end

# Lagrange interpolant from barycentric form.
# Algorithm 31.

function LagrangeInterpolant(x,xvalues::Array{ComplexF64,1},fvalues::Array{ComplexF64,1},wvalues::Array{ComplexF64,1})
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

function LagrangeInterpolatingPolynomials(x::ComplexF64,xvalues::Array{ComplexF64,1},wvalues::Array{ComplexF64,1})
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


# Interpolation from a course to a fine grid in 2D.
# Algorithm 35.


# # # INCOMPLETE # # #
#=
function CourseToFineInterpolation2D(x::Array{ComplexF64,1},y::Array{ComplexF64,1}, 
    f::Array{ComplexF64,1},ξ::Array{ComplexF64,1},η::Array{ComplexF64,1})
        

    Nₒ = length(x) - 1
    Mₒ = length(y) - 1
    @assert((Nₒ,Mₒ = size(f)))

    Nₙ = length(ξ) - 1
    Mₙ = length(η) - 1


    w = BarycentricWeights(x)
    T = PolymialInterpolationMatrix(x,w,ξ)
    F̅ = 

end
=#


#
# Algorithm 36.

function LagrangeInterpolantDerivative(x,xvec::Array{ComplexF64,1},fvec::Array{ComplexF64,1},wvec::Array{ComplexF64,1})
    @assert(length(xvec) == length(fvec) == length(wvec))
    N = length(xvec) - 1
    atNode = false
    numerator = 0
    denominator = 0
    i = Int
    for j in 0:N
        if isapprox(x,xvec[j+1]) == true
            atNode = true
            p = fvec[j+1]
            denominator = -wvec[j+1]
            i = j
        end
    end

    if atNode == true
        for j in 0:N
            if j != i
                numerator = numerator + wvec[j+1]*(p-fvec[j+1])/(x-xvec[j+1])
            end
        end

    else
        p = LagrangeInterpolation(x,xvec,fvec,wvec)
        for j in 0:N
            t = wvec[j+1]/(x-xvec[j+1])
            numerator = numerator + t*(p-fvec[j+1])/(x-xvec[j+1])
            denominator = denominator + t
        end   
    end

    return numerator/denominator
end


#
# Algorithm 37.

function PolynomialDerivateMatrix(xvec::Array{ComplexF64,1})
    N = length(xvec) - 1
    D = Array{ComplexF64,2}(fill(0,(N+1,N+1)))
    wvec = BarycentricWeights(xvec)

    for i in 0:N
        for j in 0:N
            if i != j
                D[i+1,j+1] = wvec[j+1]/(wvec[i+1]*(xvec[i+1] - xvec[j+1]))
                D[i+1,i+1] = -D[i+1,j+1]
            end
        end
    end
    
    return D
end


#
# Algorithm 38.

function mthOrderPolynomialDerivativeMatrix(m::Int, xvec::Array{ComplexF64,1})
    N = length(xvec) - 1

    wvec = BarycentricWeights(xvec)
    Dₘ = PolynomialDerivateMatrix(xvec)

    if m == 1
        return Dₘ
    end

    for k in 2:m
        for i in 0:N
            Dₘ[i+1,i+1] = 0
            for j in 0:N
                if j != i
                    Dₘ[i+1,j+1] = k/(xvec[i+1] - xvec[j+1]) * (wvec[j+1]/wvec[i+1] * Dₘ[i+1,i+1] - Dₘ[i+1,j+1])
                    Dₘ[i+1,i+1] = Dₘ[i+1,i+1] - Dₘ[i+1,j+1]
                end
                
            end
        end
    end

    return Dₘ
end


#
# Algorithm 39.



#
# Algorithm 40.

function FastChebyshevDerivative(f::Array{ComplexF64,1})
    N = length(f) - 1
    C,S = InitializeFCT(N)
    f̃ = FastChebyshevTransform(f,C,S,1)
    f̃¹ = ChebyshevDerivativeCoefficients(f̃)
    𝔇f = FastChebyshevTransform(f̃¹,C,S,-1)

    return 𝔇f
end




# # # TESTING ENVIRONMENT # # #
#=
test1 = Array{ComplexF64,1}(randn(2^5 + 1))
test2 = Array{ComplexF64,1}(randn(2^5 + 1))
weighttest1 = BarycentricWeights(test1)
T = PolymialInterpolationMatrix(test1,weighttest1,test2)
size(T)

mthOrderPolynomialDerivativeMatrix(4,test1)

C,S = InitializeFCT(2^5)
FastCosineTransform(test1,C,S,1)

FastChebyshevDerivative(test1)
=#


