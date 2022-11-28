
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


# Values.

v = 0.2

# Truncation of the exact solution for the advection-diffusion equation.
M = 63

function phi(x,t)
    s = 0
    for k in (-M):M
        s = s + 1/(2^(abs(k))) * exp(im*k*(x-t)/2 - v*k*k*t/2) * exp(im*k*x/2)
    end
    #return 1/(2^(abs(M))) * exp(im*M*(x-t)/2 -v*M*M*t) * exp(im*M*x/2)
    
    return real(s)
end


function phiInitial(x)
    return 3/(5-4*cos(x))
end


#Computing the time derivative.
function CollocationTimeDerivative(u::AbstractArray{Float64,1},D::AbstractArray{Float64,2})

    return -(1/2)*D*(u-v*(1/2)*D*u)
    #return -D*(u-v*D*u)
end

#Integrating
#=
function CollocationRK3()
    
end
=#


N = 64
D = FourierDerivativeMatrix(N)
f(u,p,t) = CollocationTimeDerivative(u,D)
u0 = ComputeValues(N,phiInitial)
tspan = (0.0,1.0)
prob = ODEProblem(f,u0,tspan)
sol = solve(prob,Tsit5())



scatter(ComputeNodes(N),sol[13],ylims = (0,3))
plot!(ComputeNodes(128),real(phi.(ComputeNodes(128),sol.t[13])))


plot!(ComputeNodes(128),real(phi.(ComputeNodes(128),sol.t[1000])))
plot!(ComputeNodes(128),real(phi.(ComputeNodes(128),sol.t[2049])))


uapprox1000(x) = FourierInterpolantFromNodes(x,ComputeNodes(32),sol[1000])

ComputeNodes(16)


xs = range(-1,1,length=100)
ts = range(0,10.0,length=1000)
surface(xs,ts,phi,xlabel = "x",ylabel = "t")



plot(xs,phiInitial.(xs))