using LinearAlgebra
using SparseArrays
using Plots
using FastGaussQuadrature
using FFTW
using DifferentialEquations
include("DiscreteFourierCoefficients.jl")
include("OtherDFTAlgorithms.jl")
include("DiscreteFourierTransform.jl")
include("LegendreAndChebyshev.jl")



k = 4.0
Re = 1.0
function TaylorGreen(t,x,y)
    return 2*k*cos(k*x)*cos(k*y)*exp(-2*k*k*t/Re)
end



function FourierInterpolantFromModes2D(fcoefficients::Array{ComplexF64,2},x::Float64,y::Float64)
    Nx,Ny = size(fcoefficients)
    NxDiv2 = Int(Nx/2)
    NyDiv2 = Int(Ny/2)

    s=0

    for m in (-NxDiv2):(NxDiv2-1)
        for n in (-NyDiv2):(NyDiv2-1)
            s += fcoefficients[m+NxDiv2+1,n+NyDiv2+1]*exp(im*(m*x+n*y))     
        end
    end

    return real(s)
end



# Assumes square. Fixed switchup of indices.
function AuxMatrices(M::Int)

    Mdiv2 = Int(M/2)
    AuxMat1 = AbstractArray{ComplexF64,2}(fill(0,M,M))
    AuxMat2 = AbstractArray{ComplexF64,2}(fill(0,M,M))
    AuxMat3 = AbstractArray{Float64,2}(fill(0,M,M))
    AuxMat4 = AbstractArray{Float64,2}(fill(0,M,M))

    for j in (-Mdiv2):(Mdiv2-1)
        for k in (-Mdiv2):(Mdiv2-1)
            AuxMat1[j+Mdiv2+1,k+Mdiv2+1] = im*k
            AuxMat2[j+Mdiv2+1,k+Mdiv2+1] = im*j

            AuxMat3[j+Mdiv2+1,k+Mdiv2+1] = -(j*j+k*k)
            AuxMat4[j+Mdiv2+1,k+Mdiv2+1] = 1/(j*j+k*k)
        end
    end

    # To avoid division by zero.
    AuxMat4[Mdiv2+1,Mdiv2+1] = 0

    return AuxMat1,AuxMat2,AuxMat3,AuxMat4
end

# Also assumes square.
function ZeroPad(w::AbstractArray{ComplexF64,2})
    N,N = size(w)
    M = Int(3/2 * N)
    
    wPadded = Array{ComplexF64,2}(fill(0,M,M))
    for m in 1:N
        for n in 1:N
            wPadded[m,n] = w[m,n]
        end
    end

    return wPadded
end


#Assumes that the auxiliary matrices are of length M = 3N/2. Also removes padding.
function TurbulenceEquationTimeDerivative(w::AbstractArray{ComplexF64,2},AuxMat1::AbstractArray{ComplexF64,2},AuxMat2::AbstractArray{ComplexF64,2},AuxMat3::AbstractArray{Float64,2},AuxMat4::AbstractArray{Float64,2},N::Int)

    w = ZeroPad(w)
    w = 4/(9*N*N)*fftshift(fft(ifft((w.*AuxMat4).*AuxMat1).*ifft(w.*AuxMat2) - ifft(w.*AuxMat1).*ifft((w.*AuxMat4).*AuxMat2))) + AuxMat3.*w
    # Removes zero-padding.
    wTimeDerivative = AbstractArray{ComplexF64,2}(fill(0,N,N))
    for m in 1:N
        for n in 1:N
            wTimeDerivative[m,n] = w[m,n]
        end
    end

    return wTimeDerivative   
end


# Uses Tsit5(). We fix M = 3N/2.
function TurbulenceEquationDriver(wInitial::AbstractArray{ComplexF64,2},N::Int)

    AuxMat1,AuxMat2,AuxMat3,AuxMat4 = AuxMatrices(Int(3*N/2))

    g(w,p,t) = TurbulenceEquationTimeDerivative(w,AuxMat1,AuxMat2,AuxMat3,AuxMat4,N)
    tspan = (0.0,1.0)
    prob = ODEProblem(g,wInitial,tspan)
    sol = solve(prob,Tsit5())

    return sol
end


# Initial condition for the modal coefficients for the Taylor-Green vortex. Not zero-padded.
function TaylorGreenInitial(N::Int)
    xs = ys = ComputeNodes(N,0,2*pi)
    InitialCondition = 1/(N*N) * (fft([TaylorGreen(0.0,x,y) for x=xs,y=ys]))

    return InitialCondition
end


# Again assumes that we want norm of a square array of values. t indicates the time
# that we are taking our true function.
function ComputeL2Error(f::Function,t::Float64,g::Function,N::Int)
    s=0

    for i in 1:N
        for j in 1:N
            s += (f(t,2*pi*i/N,2*pi*j/N) - g(2*pi*i/N,2*pi*j/N))^2
        end
    end
    
    s = 1/(N*N) * sqrt(s)

    return s
end

# Assumes again that matrices are square.
function ComputeLinfError(f::Function,t::Float64,g::Function,N::Int)

    Mat = AbstractArray{Float64,2}(fill(0,N,N))
    for i in 1:N
        for j in 1:N
            Mat[i,j] = abs( f(t,2*pi*i/N,2*pi*j/N) - g(2*pi*i/N,2*pi*j/N) )
        end
    end

    return maximum(Mat)
end



### TESTING ENVIRONMENT ###
#heatmap(A/(maximum(A)),clims=(-1,1))

xs = ys = range(0,2*pi,length=128)
A = [TaylorGreen(0.0,xs[i],ys[j]) for i in eachindex(xs),j in eachindex(ys)]
heatmap(xs,ys,A,c=cgrad([:blue,:white,:red]))


wInitial = TaylorGreenInitial(128)
w0pad = ZeroPad(fftshift(wInitial))
maximum(real(wInitial))
maximum(imag(wInitial))

AuxMat1,AuxMat2,AuxMat3,AuxMat4 = AuxMatrices(Int(3*128/2))

test = AuxMat3.*w0pad

maximum(real(test))


test2 = 4/(9*128*128)*fftshift(fft(ifft((w0pad.*AuxMat4).*AuxMat1).*ifft(w0pad.*AuxMat2) - ifft(w0pad.*AuxMat1).*ifft((w0pad.*AuxMat4).*AuxMat2)))

maximum(real(test2))
maximum(imag(test2))



fInitial(x,y) = FourierInterpolantFromModes2D(fftshift(wInitial),x,y)

B = [fInitial(xs[i],ys[j]) for i in eachindex(xs),j in eachindex(ys)]
heatmap(xs,ys,B-A,c=cgrad([:blue,:white,:red]))



test = TurbulenceEquationTimeDerivative(wInitial,t1,t2,t3,t4,128)
maximum(imag(test))
maximum(real(test))


test2 = t3.*ZeroPad(wInitial)
maximum(real(test2))

maximum(real(wInitial))


w = TaylorGreenInitial(128)
out = TurbulenceEquationDriver(w,128)


