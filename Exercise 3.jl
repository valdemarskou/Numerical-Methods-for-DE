using LinearAlgebra
using SparseArrays
using Plots
using FFTW
using DifferentialEquations
include("DiscreteFourierCoefficients.jl")
include("OtherDFTAlgorithms.jl")
include("DiscreteFourierTransform.jl")





function TaylorGreen(t,x,y,Re)
    k = 4.0
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




# Everything is assumed zero-padded.
function TurbulenceEquationTimeDerivative(w::AbstractArray{ComplexF64,2},AuxMat1::AbstractArray{ComplexF64,2},AuxMat2::AbstractArray{ComplexF64,2},AuxMat3::AbstractArray{Float64,2},AuxMat4::AbstractArray{Float64,2},N::Int)
    wTimeDeriv =(1/Re) * AuxMat3.*w
    wTimeDeriv += fft(((w.*AuxMat4).*AuxMat1).*ifft(w.*AuxMat2) - ifft(w.*AuxMat1).*ifft((w.*AuxMat4).*AuxMat2))

    return wTimeDeriv  
end

function TurbulenceEquationDriver(wInitial::AbstractArray{ComplexF64,2},t::Float64,Re::Int)
    N,N = size(wInitial)

    AuxMat1,AuxMat2,AuxMat3,AuxMat4 = AuxMatrices(N)
    wInitial = ZeroPad(wInitial)
    AuxMat1 = ZeroPad(AuxMat1)
    AuxMat2 = ZeroPad(AuxMat2)
    AuxMat3 = ZeroPad(AuxMat3)
    AuxMat4 = ZeroPad(AuxMat4)

    g(u,p,t) = TurbulenceEquationTimeDerivativeNoPad(u,AuxMat1,AuxMat2,AuxMat3,AuxMat4,Re)
    tspan = (0.0,t)
    prob = ODEProblem(g,wInitial,tspan,Re)
    sol = solve(prob,Tsit5())

    return sol
end


# Does not zero-pad.
function TurbulenceEquationDriverNoPad(wInitial::AbstractArray{ComplexF64,2},t::Float64,Re::Int)

    N,N = size(wInitial)

    AuxMat1,AuxMat2,AuxMat3,AuxMat4 = AuxMatrices(N)
    #AuxMat1 = ZeroPad(AuxMat1)
    #AuxMat2 = ZeroPad(AuxMat2)
    #AuxMat4 = ZeroPad(AuxMat4)

    g(u,p,t) = TurbulenceEquationTimeDerivativeNoPad(u,AuxMat1,AuxMat2,AuxMat3,AuxMat4,Re)
    tspan = (0.0,t)
    prob = ODEProblem(g,wInitial,tspan,Re)
    sol = solve(prob,Tsit5())

    return sol
end

function TurbulenceEquationTimeDerivativeNoPad(w::AbstractArray{ComplexF64,2},AuxMat1::AbstractArray{ComplexF64,2},AuxMat2::AbstractArray{ComplexF64,2},AuxMat3::AbstractArray{Float64,2},AuxMat4::AbstractArray{Float64,2},Re::Int)
    wTimeDeriv =(1/Re) * AuxMat3.*w
    wTimeDeriv += fft(((w.*AuxMat4).*AuxMat1).*ifft(w.*AuxMat2) - ifft(w.*AuxMat1).*ifft((w.*AuxMat4).*AuxMat2))

    return wTimeDeriv
end


# Initial condition for the modal coefficients for the Taylor-Green vortex. Not zero-padded.
function TaylorGreenInitialModes(N::Int,t::Float64)
    xs = ys = ComputeNodes(N,0,2*pi)
    InitialCondition = 1/(N*N) * (fft([TaylorGreen(t,x,y,1.0) for x=xs,y=ys]))

    return InitialCondition
end


# Again assumes that we want norm of a square array of values. t indicates the time
# that we are taking our true function.
function ErrorL2(f::Function,t::Float64,w::AbstractArray{ComplexF64,2},N::Int)
    s=0

    for i in 1:N
        for j in 1:N
            s += (f(t,2*pi*i/N,2*pi*j/N,1) - FourierInterpolantFromModes2D(w,2*pi*i/N,2*pi*j/N))^2
        end
    end
    
    s = 1/(N*N) * sqrt(s)

    return s
end

# Assumes again that matrices are square.
function ErrorLinf(f::Function,t::Float64,w::AbstractArray{ComplexF64,2},N::Int)
    w = fftshift(w)
    Mat = AbstractArray{Float64,2}(fill(0,N,N))
    for i in 1:N
        for j in 1:N
            Mat[i,j] = abs( f(t,2*pi*i/N,2*pi*j/N) - FourierInterpolantFromModes2D(w,2*pi*i/N,2*pi*j/N) )
        end
    end

    return maximum(Mat)
end


function LinfNormOfInterPolant(w::AbstractArray{ComplexF64,2})
    w = fftshift(w)

    xs = ys = range(0,2*pi,length=size(w,1))
    Mat = [abs(FourierInterpolantFromModes2D(w,x,y)) for x=xs, y = ys]

    return maximum(Mat)
end

# Computes the modal coefficients for ψ, based on the computed coefficients for w.
function StreamCoefficients(w::AbstractArray{ComplexF64,2})
    a4 = AuxMatrices(size(w,1))[4]

    return (a4.*w)
end



# Computes the modal coefficients for w, based on the randomly generated coefficients for ψ.
function FluidCoefficients(ψ::AbstractArray{ComplexF64,2})
    a3 = AuxMatrices(size(w,1))[3]
    return -(a3.*ψ)
end

function RandomFluid(N::Int)
    ψ = randn(Complex{Float64}, N,N)
    ψ = ψ * 1/(sqrt(2*TotalEnergy(ψ)))

    return ψ
end


#Take input in the nyquist form. Assumes N_x = N_y.
function EnergySpectrum(w::AbstractArray{ComplexF64,2},k::Int)
    ψ = fftshift(StreamCoefficients(w))
    Nx,Ny = size(ψ)
    NxDiv2 = Int(Nx/2)
    NyDiv2 = Int(Ny/2)

    t = 0

    for m in (-NxDiv2):(NxDiv2-1)
        for n in (-NyDiv2):(NyDiv2-1)
            if(k*k < m*m+n*n && m*m+n*n < (k+1)*(k+1))
                t += real(ψ[m+NxDiv2+1,n+NxDiv2+1])^2 + imag(ψ[m+NxDiv2+1,n+NxDiv2+1])^2
            end    
        end
    end
    
    return t
end


function TotalEnergy(w::AbstractArray{ComplexF64,2})
    return sum(real(w).*real(w)) + sum(imag(w).*imag(w))
end


function DipoleInitialCondition(N::Int)
    d = pi
    xs = ys = ComputeNodes(N,0,2*pi)
    A = [-5*exp(-((x-d)^2 + (y-21/20 *d)^2)/(0.2*d)) for x = xs, y = ys] 
    B = [5*exp(-((x-d)^2 + (y-19/20 *d)^2)/(0.2*d)) for x = xs, y = ys]
    
    return 1/(N*N) * fft(A+B)
    #return A+B
end


# For plotting purposes.
function HeatMapOfInterpolant(w::AbstractArray{ComplexF64,2},r)
    w = fftshift(w)
    xs = ys = range(0,2*pi,length = size(w,1))
    Mat = [FourierInterpolantFromModes2D(w,x,y) for x = xs, y = ys]
    heatmap(xs,ys,Mat,c=:vik,clims=(-r,r))
end

function HeatMapOfInterpolant(w::AbstractArray{ComplexF64,2})
    w = fftshift(w)
    e = LinfNormOfInterPolant(w)
    xs = ys = range(0,2*pi,length = size(w,1))
    Mat = [FourierInterpolantFromModes2D(w,x,y) for x = xs, y = ys]

    Mat = Mat * 1/(maximum(abs.(Mat)))
    heatmap(xs,ys,Mat,c=:vik)
end


### TESTING ENVIRONMENT ###

# To run the simulations, you first need to generate an initial field. 
# There is TaylorGreenInitialModes(N,0.0), which takes integer input and returns
# an NxN matrix of modal coefficients for the classical Taylor Green Vortex.
# There is also DipoleInitialCondition(N), which is just a different and fun 
# possible initial field. Finally, to generate a random field we use
# RandomFluid(N) which again creates a NxN Matrix.
# After this, use TurbulenceEquationDriver or TurbulenceEquationDriverNoPad, 
# which takes as input the initial field, the stopping time t, 
# and the reynolds Number of the flow. Then the result can be plotted using
# HeatMapOfInterpolant, which takes as input a set of moedal coefficients at 
# a desired timestep.

#Example: 

InitialField = TaylorGreenInitialModes(64,0.0)
solution = TurbulenceEquationDriverNoPad(InitialField, 0.1,1)



# If using the zero-padded driver, always use the RemovePad() function on the 
# final output. And remember to use fftshift to convert between normal and
# nyquist frequencies in the modes.



# Plots for a simulation of a random initial fluid at some different timesteps.
#=
w = RandomFluid(128)
sol = TurbulenceEquationDriver(w,20.0,1000)

p1 = HeatMapOfInterpolant(RemovePad(sol[1]))
p1 = plot(p1, xlabel = "t = 0.0")
p2 = HeatMapOfInterpolant(RemovePad(sol[20]))
p2 = plot(p2,xlabel = "t = 5.9")
p3 = HeatMapOfInterpolant(RemovePad(sol[40]))
p3 = plot(p3, xlabel = "t = 14.4")
p4 = HeatMapOfInterpolant(RemovePad(sol[54]))
p4 = plot(p4, xlabel = "t = 20.0")
out = plot(p1,p2,p3,p4,layout = (2,2))
#savefig(out,"RandomFluid.png")
=#



# Plot that measures the computation time. Will take a few minutes to run.
#=
ComputationTimes = Array{Float64,1}(fill(0,33))
L2Errors = Array{Float64,1}(fill(0,33))
LinfErrors = Array{Float64,1}(fill(0,33))
for k in 32:64
    t = TaylorGreenInitialModes(2*k,0.0)
    sol = @timed TurbulenceEquationDriver(t,0.1)

    ComputationTimes[k-32+1] = sol[2]
    L2Errors[k-32+1] = ComputeL2Error(TaylorGreen,0.1,fftshift(last(sol[1])),2*k)
    LinfErrors[k-32+1] = ComputeLinfError(TaylorGreen,0.1,fftshift(last(sol[1])),2*k)
end
=#
#=
evens = AbstractArray{Float64,1}(fill(0,33))
for k in 0:32
    evens[k+1] = 64+2*k
end

p1 = plot(evens,ComputationTimes,linestyle=:dash,title="Computation Time",xlabel="N",ylabel = "time (s)",label="Computation Time",color =:red)
p2 = scatter(evens,L2Errors,marker = :cross,color=:darkblue,yaxis = :log, label = "L2-Error",xlabel = "N",title = "L2 and L∞ Errors",ylabel = "Error")
scatter!(evens,LinfErrors,marker = :cross, color=:turquoise,yaxis =:log,label = "L∞-Error",xlabel = "N")
fig = plot(p1,p2,layout=(1,2))
savefig(fig,"ComputationTimeAndErrors.png")
=#



# Some other plot.
#=
w = TaylorGreenInitialModes(128,0.0)
sol = TurbulenceEquationDriver(w,0.1,1)

E1=[EnergySpectrum(sol[1],k) for k in 1:135]
E2 = [EnergySpectrum(last(sol),k) for k in 1:135]
asym = [(1e-25)/(k^3) for k in 1:135]



p1 = scatter(1:100,E1,yaxis=:log,xaxis=:log,label ="t = 0",marker=:cross,color=:darkred,title="Energy Spectrum For Different time")
scatter!(1:100,E2,label ="t = 0.1",marker=:cross,color=:pink)

plot!(asym,linestyle=:dash,color=:black,label="1/k³")
#savefig(p1, "spectrumplot.png")
=#