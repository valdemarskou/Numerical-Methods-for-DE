include("LegendreAndChebyshev.jl")


#The Fourier Collocation Time Derivative for the Advection-Diffusion Equation. Here
# ν>0 is assumed to be constant. Algorithm 41.

function FourierCollocationTimeDerivative(Φ::Array{Float64,1},D::Array{Float64,2},v::Float64)

    N = length(Φ)
    F = MxVDerivative(D,Φ,0)

    for j in 0:(N-1)
        F[j+1] = v*F[j+1] - Φ[j+1]
    end

    PhiDeriv = MxVDerivative(D,F,0)

    return PhiDeriv
end


# Fourier Collocation Time Derivative for the Advection-Diffusion equation, using the
# fast fourier transform instead of the matrix Multiplication trick.
# exercise 4.1 (alternate algorithm 41)

function FourierCollocationTimeDerivative(Φ::Array{Float64,1},v::Float64)
 
    N = length(Φ)
    DΦ = FourierDerivativeByFFT(Φ,1)
    


    
end


#
# Algorithm 42.

# # # INCOMPLETE DUE TO NOT KNOWING HOW TO DEFINE THE G ARRAY # # #

#=
function CollocationStepByRK3(t_n::Float64, Δt::Float64,Φ::Array{Float64,1},D::Array{Float64,2},v::Float64)
    a = [0, -5/9, -153/128]
    b = [0, 1/3, 3/4]
    g = [1/3, 15/16, 8/15]
    N = length(Φ)

    for m in 1:3
        t = t_n + b[m]*Δt
        PhiDeriv = FourierCollocationTimeDerivative(Φ,D,v)
        for j in 0:(N-1)
            
        end
    end
end
=#


# Advection-Diffusion Time Derivative for Fourier Galerkin. Requires N to be even.
# Algorithm 44.

function AdvectionDiffusionTimeDerivative(PhiCoeffs::Array{Float64,1},v::Float64)
    N = length(PhiCoeffs) - 1
    Ndiv2 = Int(N/2)
    PhiDerivCoeffs = Array{ComplexF64,1}(fill(0,N+1))
    
    for k in (-Ndiv2):(Ndiv2)
        PhiDerivCoeffs[k+1+Ndiv2] = -(im*k + v*k^2)*PhiCoeffs[k+1+Ndiv2]
    end

    return PhiDerivCoeffs
end


# Take One Time Step of the Fourier Galerkin Method. For both real and complex array.
# Algorithm 45.

function FourierGalerkinStep(t_n::Float64, Δt::Float64,PhiCoeffs::Array{Float64,1},v::Float64)

    a = [0, -5/9, -153/128]
    b = [0, 1/3, 3/4]
    g = [1/3, 15/16, 8/15]

    N = length(PhiCoeffs) - 1
    Ndiv2 = N/2

    for m in 1:3
        t = t_n + Δt * b[m]
        PhiDerivCoeffs = AdvectionDiffusionTimeDerivative(PhiCoeffs,v)
        G = 0
        for k in (-Ndiv2):(Ndiv2)
            G = a[m]*G + PhiDerivCoeffs[k+1+Ndiv2]
            PhiCoeffs[k+1+Ndiv2] = PhiCoeffs[k+1+Ndiv2] + g[m]*(Δt)*G 
        end
    end

    return PhiCoeffs
end

function FourierGalerkinStep(t_n::Float64, Δt::Float64,PhiCoeffs::Array{ComplexF64,1},v::Float64)

    a = [0, -5/9, -153/128]
    b = [0, 1/3, 3/4]
    g = [1/3, 15/16, 8/15]

    N = length(PhiCoeffs) - 1
    Ndiv2 = N/2

    for m in 1:3
        t = t_n + Δt * b[m]
        PhiDerivCoeffs = AdvectionDiffusionTimeDerivative(PhiCoeffs,v)
        G = 0
        for k in (-Ndiv2):(Ndiv2)
            G = a[m]*G + PhiDerivCoeffs[k+1+Ndiv2]
            PhiCoeffs[k+1+Ndiv2] = PhiCoeffs[k+1+Ndiv2] + g[m]*(Δt)*G 
        end
    end

    return PhiCoeffs
end



# Direct Synthesis of the Fourier Galerkin Solution.
# Algorithm 46.

function EvaluateFourierGalerkinSolution(x,PhiCoeffs)

    N = length(PhiCoeffs) + 1
    Ndiv2 = Int(N/2)
    Phi = 0

    for k in (-Ndiv2):(Ndiv2)
        Phi = Phi + PhiCoeffs[k+1+N/2]*exp(im*k*x)
    end

    return Phi
end


# A Driver for the Fourier Galerkin Approximation.
# Algoritm 47.

function FourierGalerkinDriver(N_T::Int,T::Float64,N_out::Int,PhiInitialCoeffs::Array{ComplexF64,1},v::Float64)

    N = length(PhiInitialCoeffs) - 1
    Δt = T/(N_T)
    t_n = 0
    PhiCoeffs = PhiInitialCoeffs
    
    Phi = Array{ComplexF64,1}(fill(0,N_out + 1))

    for n in 0:(N-1)
        PhiCoeffs = FourierGalerkinStep(t_n, Δt, PhiCoeffs,v)
        t_n = (n+1)* Δt
    end

    Δx = 2*pi/(N_out)

    for j in 0:N_out
        x_j = j* Δx
        Phi[j+1] = EvaluateFourierGalerkinSolution(x_j,PhiCoeffs)
    end

    return Phi
end


# Direct (Slow) Computation of the Convolution Sum
# Algorithm 48.

function DirectConvolutionSum(V̂::Array{ComplexF64,1},Ŵ::Array{ComplexF64,1})

    @assert(length(V̂) == length(Ŵ))
    N = length(V̂) - 1
    Ndiv2 = Int(N/2)

    VWconvolution = Array{ComplexF64,1}(fill(0,N+1))

    for k in -Ndiv2:Ndiv2
        for p in max(-Ndiv2, k - Ndiv2):min(Ndiv2, Ndiv2 + k)
            VWconvolution[k+Ndiv2+1] = VWconvolution[k+Ndiv2+1] + V̂[k - p + Ndiv2 + 1]*Ŵ[p + Ndiv2 + 1]
        end
        
    end

    return WVconvolution
end


# Computation of the Convolution Sum with the FFT.
# Algorithm 49.

function FastConvolutionSum(V̂::Array{ComplexF64,1},Ŵ::Array{ComplexF64,1})

    @assert(length(V̂) == length(Ŵ))
    N = length(V̂) - 1
    Ndiv2 = Int(N/2)
    M = 2*N

    Ṽ = Array{ComplexF64,1}(fill(0,M))
    W̃ = Array{ComplexF64,1}(fill(0,M))

    for k in 0:Ndiv2
        Ṽ[k + 1] = V̂[k + Ndiv2 + 1]
        W̃[k + 1] = W̃[k + Ndiv2 + 1]
    end

    for k in reverse(-Ndiv2:(-1))
        Ṽ[M + k + 1] = V̂[k + Ndiv2 + 1]
        W̃[M + k + 1] = V̂[k + Ndiv2 + 1]
    end

    V = bfft(Ṽ)
    W = bfft(W̃)

    Q = Array{ComplexF64,1}(fill(0,M))

    for k in 0:(M-1)
        Q[k+1] = V[k+1]*W[k+1]
    end

    Q̃ = fft(Q)

    VWconvolution = Array{ComplexF64,1}(fill(0,N+1))

    for k in 0:Ndiv2
        VWconvolution[k+Ndiv2+1] = Q[k+1]
    end
    
    for k in reverse(-Ndiv2:(-1))
        VWconvolution[k+Ndiv2+1] = Q[M+k+1]
    end

    return VWconvolution
end

#=
test1 = Array{ComplexF64,1}(randn(2^6 + 1))
test2 = Array{ComplexF64,1}(randn(2^6 + 1))

FastConvolutionSum(test1,test2)
=#


#
# Algorithm 50.

function ()
    
end