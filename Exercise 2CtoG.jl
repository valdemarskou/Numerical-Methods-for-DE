
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

# Problem c)

L = 50
T = 1
c = 1
x0 = 0

function uInitial(x)
    return c/2 * sech(sqrt(c)/2 *(x-x0))^2
end

function uInitial2(x)
    return uInitial(L*(x/pi - 1))
end

function uExact(x,t)
    return uInitial(x-c*t)
end

function uExact2(x,t)
    return uExact(L*(x/pi - 1),t*T)
end


function KDVEquationGalerkinTimeDerivative(û::AbstractArray{ComplexF64,1},AuxVec::AbstractArray{Float64,1})


    #ûtimederiv = AbstractArray{ComplexF64,1}(fill(0,N+1))
    ûtimederiv = fftshift(fft(ifft(û) .* ifft(û.*AuxVec)))

    ûtimederiv = -6*T*(pi/L)*ûtimederiv + T*((pi/L)^3)*im*(AuxVec.*AuxVec.*AuxVec).*û

    return ûtimederiv
end


function KDVEquationGalerkinDriver(N::Int,c,x0)
    AuxVec = AbstractArray{Float64,1}(fill(0,N))
    for k in 0:(N-1)
        AuxVec[k+1] = k-N/2
    end

    g(u,p,t) = KDVEquationGalerkinTimeDerivative(u,AuxVec)
    u0 = fft(ComputeValues(N,uInitial2,0,2*pi))
    tspan = (0.0,1.0)
    prob = ODEProblem(g,u0,tspan)
    sol = solve(prob)

    return sol
end

ComputeValues(16,uInitial2,0,2*pi)



function KDVEquationCollocationTimeDerivative(u::AbstractArray{Float64,1},D::AbstractArray{Float64,2})

    #N = length(u)
    #D = FourierDerivativeMatrix(N)

    return (-6*T*(pi/L)*(u .*(D*u)) - T*((pi/(L))^(3))*D*D*D*u)

    #return (-6*T*(2*pi/L)*(u .*(FourierDerivativeByFFT(u,1))) - T*((2*pi/(L))^3)*FourierDerivativeByFFT(FourierDerivativeByFFT(FourierDerivativeByFFT(u,1),1),1))
    
end



function KDVEquationCollocationDriver(N::Int,c,x0)
    D = FourierDerivativeMatrix(N)
    u0 = uInitial2.(ComputeNodes(N),0.5,-20)
    g(u,p,t) = KDVEquationCollocationTimeDerivative(u,D)
    #tspan=(0.0,1.0)
    tspan = (0.0,1.0)
    prob = ODEProblem(g,u0,tspan)
    sol = solve(prob)
    
    return sol
end



function ComputeL2Error(N::Int,c,x0)
    out = KDVEquationCollocationDriver(N,c,x0)
    k=1
    L2Error=AbstractArray{Float64,1}(fill(0,length(out.t)))
    for t in out.t
        s=2*pi/N * sum((uExact2.(ComputeNodes(N),t,c,x0)-out[k]).*(uExact2.(ComputeNodes(N),t,c,x0)-out[k])) 
        L2Error[k] = s+0.00000000000000000000000000001
        k=k+1
    end

    return L2Error
    #plot(out.t,l2error)    
end

function ComputeLinfError(N::Int,c,x0)
    out = KDVEquationCollocationDriver(N,c,x0)
    k=1
    LinfError=AbstractArray{Float64,1}(fill(0,length(out.t)))
    for t in out.t 
        s = maximum(abs.(uExact2.(ComputeNodes(N),t,c,x0)-out[k]))
        LinfError[k] = s+0.00000000000000000000000000001
        k=k+1
    end
    return LinfError
end





#Plots
sol16 = KDVEquationCollocationDriver(16,0.5,-20)
sol64 = KDVEquationCollocationDriver(64,0.5,-20)
sol128 = KDVEquationCollocationDriver(128,0.5,-20)

l2error16 = ComputeL2Error(16,0.5,-20)
l2error64 = ComputeL2Error(64,0.5,-20)
l2error128 = ComputeL2Error(128,0.5,-20)

outa=plot(sol16.t,l2error16,yaxis=:log,linestyle=:dash,title="L2-error, c=0.5",label ="N=16")
plot!(sol64.t,l2error64,yaxis=:log,linestyle=:dash,label="N=64")
plot!(sol128.t,l2error128,yaxis=:log,linestyle=:dash,label="N=128")

l2error16b = ComputeL2Error(16,1.0,-20)
l2error64b = ComputeL2Error(64,1.0,-20)
l2error128b = ComputeL2Error(128,1.0,-20)

outb=plot(sol16.t,l2error16b,yaxis=:log,linestyle=:dash,title="L2-error, c=1.0",label ="N=16")
plot!(sol64.t,l2error64b,yaxis=:log,linestyle=:dash,label="N=64")
plot!(sol128.t,l2error128b,yaxis=:log,linestyle=:dash,label="N=128")

out=plot(outa,outb,layout=(2,1))
#savefig(out,"kdvl2errorplot.png")


#more Plots
linferror16a = ComputeLinfError(16,0.5,-20)
linferror64a = ComputeLinfError(64,0.5,-20)
linferror128a =ComputeLinfError(128,0.5,-20)

plota=plot(sol16.t,linferror16a,yaxis=:log,linestyle=:dash,title="L∞-error,c=0.5",label="N=16")
plot!(sol64.t,linferror64a,yaxis=:log,linestyle=:dash,label="N=64")
plot!(sol128.t,linferror128a,yaxis=:log,linestyle=:dash,label="N=128")


linferror16b = ComputeLinfError(16,1.0,-20)
linferror64b = ComputeLinfError(64,1.0,-20)
linferror128b = ComputeLinfError(128,1.0,-20)

plotb=plot(sol16.t,linferror16b,yaxis=:log,linestyle=:dash,title="L∞-error,c=1.0",label="N=16")
plot!(sol64.t,linferror64b,yaxis=:log,linestyle=:dash,label="N=64")
plot!(sol128.t,linferror128b,yaxis=:log,linestyle=:dash,label="N=128")


out=plot(plota,plotb,layout=(2,1))
#savefig(out,"kdvlinferror.png")


#Collision
T=120
wave1 = KDVEquationCollocationDriver(64,0.25,40)
wave2 = KDVEquationCollocationDriver(64,0.5,15)

scatter(ComputeNodes(64),wave1[1]+wave2[1],xlims=(0,2*pi))

#Computationtimeplot



ComputationTime=AbstractArray{Float64,1}(fill(0,57))

for k in 0:56
    t = @timed KDVEquationCollocationDriver(16+2*k,0.5,-20)
    ComputationTime[k+1] = t.time
end

ComputationTime
evens = AbstractArray{Float64,1}(fill(0,57))
for k in 0:56
    evens[k+1] = 16+2*k
end

out=plot(evens,ComputationTime,linestyle=:dash,title="Computation Time",xlabel="N",label="Computation Time (s)")
#savefig(out,"computationtime.png")





## TESTING ENVIRONMENT ##





M = 64

output = KDVEquationGalerkinDriver(M,1,0)



p1 = plot(ComputeNodes(512),uExact2.(ComputeNodes(512),0))
interpt0(x) = FourierInterpolantFromModes(x,output[1])
scatter!(ComputeNodes(512),interpt0.(ComputeNodes(512)))

out = plot(xs,interpt0.(xs))