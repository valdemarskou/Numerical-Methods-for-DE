
using LinearAlgebra
using SparseArrays
using Plots
using FastGaussQuadrature
using FFTW
include("DiscreteFourierCoefficients.jl")
include("OtherDFTAlgorithms.jl")
include("DiscreteFourierTransform.jl")
include("LegendreAndChebyshev.jl")



function v(x,ϵ)
    return (exp(-x/ϵ)+(x-1)-exp(-1/ϵ)*x)/(exp(-1/ϵ)-1)
end

function u(x,ϵ)
    return v((x+1)/2,ϵ)
end

# Problem a) Implementation of the LTM sparse matrix.

function LTMSolver(N::Int,ϵ::Float64) 
    A = AbstractArray{Float64,2}(fill(0,(N+1,N+1)))
    for n in 1:(N+1)
        A[1,n] = 1
        A[2,n] = (-1)^(n)
    end
    
    for n in 2:(N-2)
        A[n+1,n+1] = 1
        A[n+1,n] = 1/(2*ϵ*(2*n-1))
        A[n+1,n+2] = - 1/(2*ϵ*(2*n+3)) 
    end
    
    A[N+1,N] = 2*(2*N-3)
    A[N+1,N+1] = 4*ϵ*(2*N-3)*(N-1)
    A[N,N+1] = 2*(2*N-1)
    
    
    f̂ = AbstractArray{Float64,1}(fill(0,(N+1)))
    f̂[1] = -1
    
    ĝ = AbstractArray{Float64,1}(fill(0,(N+1)))
    
    for n in 2:(N-2)
        ĝ[n+1] = 1/(4*ϵ*(2*n-1)) * (f̂[n-1]/(2*n-3) - f̂[n+1]/(2*n+1)) - 1/(4*ϵ*(2*n+3)) * (f̂[n+1]/(2*n+1) - f̂[n+3]/(2*n+5))
    end
    
    ĝ[N] = f̂[N]
    ĝ[N+1] = f̂[N-1]

    return A\ĝ
end

# Problem a) Legendre approximation based on coefficient output.

function uapprox(x,û)
    N = length(û) - 1
    s = 0
    for n in 0:N
        s = s+û[n+1]*LegendrePolynomial(n,x)
    end
    return s
end


# Problem a) LTM final plot.
let 
    ϵ = 0.1
    û16 = LTMSolver(16,ϵ)
    û32 = LTMSolver(32,ϵ)
    û128 = LTMSolver(128,ϵ)
    
    uapprox16e(x) = abs(u(x,ϵ) - uapprox(x,û16))+0.0000000000000000001
    uapprox32e(x) = abs(u(x,ϵ) - uapprox(x,û32))+0.0000000000000000001
    uapprox128e(x) = abs(u(x,ϵ) - uapprox(x,û128))+0.0000000000000000001
    
    #almostzero = AbstractArray{Float64,1}(fill(0.0000000000000000001,100))
    
    xs = range(-1,1,length = 100)
    
    p1 = plot(xs,u.(xs,ϵ),label = "ϵ = 0.1",color = :black,ylims = (0,1))

    #p2 = plot(xs,uapprox16e.(xs),yaxis = :log, label = "N = 16",color = :blue)
    p2 = plot(xs,uapprox16e.(xs),yaxis = :log, label = "N = 16",color = :red)
    plot!(xs,uapprox32e.(xs),yaxis = :log, label = "N = 32",color = :blue)
    plot!(xs,uapprox128e,yaxis = :log,label = "N = 128",color = :green)
    
    
    ϵ = 0.01
    û16 = LTMSolver(16,ϵ)
    û32 = LTMSolver(32,ϵ)
    û128 = LTMSolver(128,ϵ)
    
    uapprox16e(x) = abs(u(x,ϵ) - uapprox(x,û16))+0.0000000000000000001
    uapprox32e(x) = abs(u(x,ϵ) - uapprox(x,û32))+0.0000000000000000001
    uapprox128e(x) = abs(u(x,ϵ) - uapprox(x,û128))+0.0000000000000000001
    
    p3 = plot(xs,u.(xs,ϵ),label = "ϵ = 0.01",title = "True solution v(x,ϵ)",color = :black,ylims = (0,1))
    p4 = plot(xs,uapprox16e.(xs),yaxis = :log, label = "N = 16",color = :red,title = "Error")
    plot!(xs,uapprox32e.(xs),yaxis = :log, label = "N = 32",color = :blue)
    plot!(xs,uapprox128e,yaxis = :log,label = "N = 128",color = :green)
    
    
    ϵ = 0.001
    û16 = LTMSolver(16,ϵ)
    û32 = LTMSolver(32,ϵ)
    û128 = LTMSolver(128,ϵ)
    
    uapprox16e(x) = abs(u(x,ϵ) - uapprox(x,û16))+0.0000000000000000001
    uapprox32e(x) = abs(u(x,ϵ) - uapprox(x,û32))+0.0000000000000000001
    uapprox128e(x) = abs(u(x,ϵ) - uapprox(x,û128))+0.0000000000000000001
    
    p5 = plot(xs,u.(xs,ϵ),label = "ϵ = 0.001",color = :black,ylims = (0,1))
    p6 = plot(xs,uapprox16e.(xs),yaxis = :log, label = "N = 16",color = :red)
    plot!(xs,uapprox32e.(xs),yaxis = :log, label = "N = 32",color = :blue)
    plot!(xs,uapprox128e,yaxis = :log,label = "N = 128",color = :green)
    
    output = plot(p1,p3,p5,p2,p4,p6,layout = (2,3))
    #savefig(output, "LTMplot.png")
    
end


function LCMSolver(N::Int,ϵ::Float64)

    D = DerivativeMatrix(0.0,0.0,gausslobatto(N+1)[1])
    A = 4*ϵ*D*D + 2*D
    for n in 1:(N+1)
        A[1,n] = 1
        A[N+1,n] = (-1)^(n)
    end

    f = AbstractArray{Float64,1}(fill(-1,N+1))
    f[1] = 0
    f[N+1] = 0
    sol = A\f
    t = sol[1]
    for n in 1:(N+1)
        sol[n] = sol[n] - t
    end

    return sol
end




# Problem a) LCM final plot.
let 
    ϵ = 0.1
    uapprox100 = LCMSolver(99,ϵ)
    uapprox200 = LCMSolver(199,ϵ)
    uapprox300 = LCMSolver(299,ϵ)
    uapprox1000 = LCMSolver(999,ϵ)

    #almostzero100 = AbstractArray{Float64,1}(fill(0.0000000000001,30))
    #almostzero200 = AbstractArray{Float64,1}(fill(0.0000000000001,60))
    #almostzero300 = AbstractArray{Float64,1}(fill(0.0000000000001,90))


    p1 = plot(gausslobatto(100),abs.(uapprox100 - u.(gausslobatto(100)[1],ϵ)), label = "N = 100",title = "ϵ = 0.1", color = :red, linestyle = :dash,xlabel = "x",ylabel = "|u(x)-Iₙu(x)|")
    plot!(gausslobatto(200),abs.(uapprox200 - u.(gausslobatto(200)[1],ϵ)),label = "N = 200", color = :blue, linestyle = :dash)
    plot!(gausslobatto(300),abs.(uapprox300 - u.(gausslobatto(300)[1],ϵ)), label = "N = 300", color = :green, linestyle = :dash)
    plot!(gausslobatto(1000),abs.(uapprox1000 - u.(gausslobatto(1000)[1],ϵ)), label = "N = 1000", color = :brown, linestyle = :dash)


    ϵ = 0.01
    uapprox100 = LCMSolver(99,ϵ)
    uapprox200 = LCMSolver(199,ϵ)
    uapprox300 = LCMSolver(299,ϵ)
    uapprox1000 = LCMSolver(999,ϵ)

    p2 = plot(gausslobatto(100),abs.(uapprox100 - u.(gausslobatto(100)[1],ϵ)), label = "N = 100",title = "ϵ = 0.01",color = :red, linestyle = :dash,xlabel = "x",ylabel = "|u(x)-Iₙu(x)|")
    plot!(gausslobatto(200),abs.(uapprox200 - u.(gausslobatto(200)[1],ϵ)),label = "N = 200", color = :blue, linestyle = :dash)
    plot!(gausslobatto(300),abs.(uapprox300 - u.(gausslobatto(300)[1],ϵ)), label = "N = 300", color = :green, linestyle = :dash)
    plot!(gausslobatto(1000),abs.(uapprox1000 - u.(gausslobatto(1000)[1],ϵ)), label = "N = 1000", color = :brown, linestyle = :dash)


    ϵ = 0.001
    uapprox100 = LCMSolver(99,ϵ)
    uapprox200 = LCMSolver(199,ϵ)
    uapprox300 = LCMSolver(299,ϵ)
    uapprox1000 = LCMSolver(999,ϵ)

    p3 = plot(gausslobatto(100),abs.(uapprox100 - u.(gausslobatto(100)[1],ϵ)), label = "N = 100",title = "ϵ = 0.001",color = :red, linestyle = :dash,xlabel = "x",ylabel = "|u(x)-Iₙu(x)|")
    plot!(gausslobatto(200),abs.(uapprox200 - u.(gausslobatto(200)[1],ϵ)),label = "N = 200", color = :blue, linestyle = :dash)
    plot!(gausslobatto(300),abs.(uapprox300 - u.(gausslobatto(300)[1],ϵ)), label = "N = 300", color = :green, linestyle = :dash)
    plot!(gausslobatto(1000),abs.(uapprox1000 - u.(gausslobatto(1000)[1],ϵ)), label = "N = 1000", color = :brown, linestyle = :dash)

    output = plot(p1,p2,p3,layout = (1,3))
    savefig(output,"LCMplot.png")
end


