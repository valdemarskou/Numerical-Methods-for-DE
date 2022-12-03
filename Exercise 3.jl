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

    ones = AbstractArray{Float64,2}(fill(1,M,M))

    AuxMat1 = im*ones.*fftfreq(M,M)
    AuxMat2 = AbstractArray{ComplexF64,2}(transpose(AuxMat1))
    AuxMat3 = -ones.*fftfreq(M,M).*(ones.*fftfreq(M,M)) - (transpose(fftfreq(M,M).*ones)).*(transpose(fftfreq(M,M).*ones))
    AuxMat4 = map((x)->-1/x,AuxMat3)

    # To avoid division by zero.
    AuxMat4[1,1] = 0.0


    return AuxMat1,AuxMat2,AuxMat3,AuxMat4
end


# Also assumes square.
function ZeroPad(w::AbstractArray{Float64,2})

    w=fftshift(w)
    N,N = size(w)
    M = Int(3/2 * N)
    
    wPadded = Array{Float64,2}(fill(0,M,M))
    for m in 1:N
        for n in 1:N
            wPadded[m,n] = w[m,n]
        end
    end

    return fftshift(wPadded)
end
function ZeroPad(w::AbstractArray{ComplexF64,2})

    w=fftshift(w)
    N,N = size(w)
    M = Int(3/2 * N)
    
    wPadded = Array{ComplexF64,2}(fill(0,M,M))
    for m in 1:N
        for n in 1:N
            wPadded[m,n] = w[m,n]
        end
    end

    return fftshift(wPadded)
end

function RemovePad(w::AbstractArray{Float64,2})

    w = fftshift(w)
    M,M = size(w)
    N = Int(2/3 * M)

    wRemoved = Array{Float64,2}(fill(0,N,N))
    for m in 1:N
        for n in 1:N
            wRemoved[m,n] = w[m,n]
        end 
    end

    return fftshift(wRemoved)
end
function RemovePad(w::AbstractArray{ComplexF64,2})

    w = fftshift(w)
    M,M = size(w)
    N = Int(2/3 * M)

    wRemoved = Array{ComplexF64,2}(fill(0,N,N))
    for m in 1:N
        for n in 1:N
            wRemoved[m,n] = w[m,n]
        end 
    end

    return fftshift(wRemoved)
end




# Removes padding. AuxMat1, AuxMat2, AuxMat4 are assumed to be zero-padded.
function TurbulenceEquationTimeDerivative(w::AbstractArray{ComplexF64,2},AuxMat1::AbstractArray{ComplexF64,2},AuxMat2::AbstractArray{ComplexF64,2},AuxMat3::AbstractArray{Float64,2},AuxMat4::AbstractArray{Float64,2},N::Int)

    wTimeDerivative = AuxMat3.*w
    w = ZeroPad(w)
    w = 4/(9*N*N)* fft(((w.*AuxMat4).*AuxMat2).*ifft(w.*AuxMat1) - ifft(w.*AuxMat2).*ifft((w.*AuxMat4).*AuxMat1))
    
    # Removes zero-padding, adds convolution terms.
    wTimeDerivative += RemovePad(w)

    return wTimeDerivative   
end


#  Stuck on evaluating...
function TurbulenceEquationDriver(wInitial::AbstractArray{ComplexF64,2})

    N,N = size(wInitial)

    AuxMat1,AuxMat2,AuxMat3,AuxMat4 = AuxMatrices(N)
    AuxMat1 = ZeroPad(AuxMat1)
    AuxMat2 = ZeroPad(AuxMat2)
    AuxMat4 = ZeroPad(AuxMat4)

    g(u,p,t) = TurbulenceEquationTimeDerivative(u,AuxMat1,AuxMat2,AuxMat3,AuxMat4,N)
    tspan = (0.0,0.1)
    prob = ODEProblem(g,wInitial,tspan)
    sol = solve(prob)

    return sol
end

function TurbulenceEquationTimeDerivativeNoPad(w::AbstractArray{ComplexF64,2},AuxMat1::AbstractArray{ComplexF64,2},AuxMat2::AbstractArray{ComplexF64,2},AuxMat3::AbstractArray{Float64,2},AuxMat4::AbstractArray{Float64,2})
    wTimeDeriv = AuxMat3.*w
    wTimeDeriv += fft(((w.*AuxMat4).*AuxMat2).*ifft(w.*AuxMat1) - ifft(w.*AuxMat2).*ifft((w.*AuxMat4).*AuxMat1))
end


# Initial condition for the modal coefficients for the Taylor-Green vortex. Not zero-padded.
function TaylorGreenInitialModes(N::Int,t::Float64)
    xs = ys = ComputeNodes(N,0,2*pi)
    InitialCondition = 1/(N*N) * (fft([TaylorGreen(t,x,y) for x=xs,y=ys]))

    return InitialCondition
end


# Again assumes that we want norm of a square array of values. t indicates the time
# that we are taking our true function.
function ComputeL2Error(f::Function,t::Float64,w::AbstractArray{ComplexF64,2},N::Int)
    s=0

    for i in 1:N
        for j in 1:N
            s += (f(t,2*pi*i/N,2*pi*j/N) - FourierInterpolantFromModes2D(w,2*pi*i/N,2*pi*j/N))^2
        end
    end
    
    s = 1/(N*N) * sqrt(s)

    return s
end

# Assumes again that matrices are square.
function ComputeLinfError(f::Function,t::Float64,w::AbstractArray{ComplexF64,2},N::Int)

    Mat = AbstractArray{Float64,2}(fill(0,N,N))
    for i in 1:N
        for j in 1:N
            Mat[i,j] = abs( f(t,2*pi*i/N,2*pi*j/N) - FourierInterpolantFromModes2D(w,2*pi*i/N,2*pi*j/N) )
        end
    end

    return maximum(Mat)
end



### TESTING ENVIRONMENT ###


xs = ys = range(0,2*pi,length=128)
A = [TaylorGreen(0.0,xs[i],ys[j]) for i in eachindex(xs),j in eachindex(ys)]

heatmap(xs,ys,A,c=:vik,clims=(-10,10))

B = [fInitial(xs[i],ys[j]) for i in eachindex(xs),j in eachindex(ys)]
heatmap(xs,ys,B-A,c=cgrad([:blue,:white,:red]))


minimum(real(TurbulenceEquationTimeDerivativeNoPad(t,a1,a2,a3,a4)))




t = TaylorGreenInitialModes(32,0.0)

ComputeL2Error(TaylorGreen,0.0,t,16)
ComputeLinfError(TaylorGreen,0.0,t,16)

maximum(real(t))
minimum(real(t))


a1,a2,a3,a4 = AuxMatrices(32)

view(fftshift(a3.*t),13,:)

view(a3,5,:)
view(t,13,:)

view(fftshift(t),5,:)


view(t.*a3,5,:)
view(t.*a3,13,:)



