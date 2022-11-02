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

v = 0.2
c = 4.0
t0 = 1.0
# Truncation of exact solution. M = ...

M = 63

function ub(x,t)
    numerator = 0.0
    denominator = 0.0

    for n in (-M:M)
        numerator += exp(-(x-2*pi*n)^2 /(4*v*t))*(x-2*pi*n)
        denominator += exp(-(x-2*pi*n)^2 /(4*v*t))
    end
    
    return (-1/(2*v*t)) * numerator/denominator
end

function uExact(x,t)
    return c + ub(x-c*t,t0+t)
end

xs = range(0,2*pi,length=100)
plot(xs,uExact.(xs,0))
plot!(xs,uExact.(xs,0.5))
plot!(xs,uExact.(xs,1))


# Computing the time derivative

function BurgersCollocationTimeDerivative(u::AbstractArray{Float64,1},D::AbstractArray{Float64,2})
    return v*D*D*u - u.*(D*u)
end




N = 64
D = FourierDerivativeMatrix(N)
f(u,p,t) = BurgersCollocationTimeDerivative(u,D)
u0 = uExact.(ComputeNodes(N),0)
tspan = (0.0,1.0)
prob = ODEProblem(f,u0,tspan)
sol = solve(prob,Tsit5())



xs = range(0,2*pi,length=100)
plot(xs,uExact.(xs,0))
#plot!(xs,uExact.(xs,0.04))
plot!(xs,uExact.(xs,sol.t[10]))
scatter!(ComputeNodes(N),sol[1])
scatter!(ComputeNodes(N),sol[10])

plot!(ComputeNodes(128),uExact.(ComputeNodes(128),sol.t[91]))



