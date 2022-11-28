
using Plots
using FastGaussQuadrature
include("DiscreteFourierCoefficients.jl")
include("OtherDFTAlgorithms.jl")
include("DiscreteFourierTransform.jl")
include("LegendreAndChebyshev.jl")


#=
function u(x)
    return 1/(2-cos(pi*x))
end
=#

# Problem a) Plot of the continuous truncation error.
let 
    N = 13
    τ = 2/(3*sqrt(3))
    ErrorVec = Array{Float64,1}(zeros(N))
    τ = τ - 1/3
    for j in 1:N
         τ = τ - 2/3 *((2-sqrt(3)))^(2*j)
        ErrorVec[j] = τ
    end
    ErrorVec
    plot2 = scatter(ErrorVec,yaxis = :log,
    xticks =[1,2,3,4,5,6,7,8,9,10,11,12,13,14],
    marker = :diamond,
    color = :blue,
    xlabel = "n",label = "τ",title = "Truncation error" )
    #savefig(plot2, "plog2.png")
end


# problem b) plot 1
let 
    function v(x)::Float64
        return 1/(2-cos(x))
    end

    function interpolant4(x)
        FourierInterpolantFromNodes(x,ComputeNodes(4),ComputeValues(4,v))
    end
    function interpolant8(x)
        FourierInterpolantFromNodes(x,ComputeNodes(8),ComputeValues(8,v))
    end
    function interpolant16(x)
        FourierInterpolantFromNodes(x,ComputeNodes(16),ComputeValues(16,v))
    end
    function interpolant32(x)
        FourierInterpolantFromNodes(x,ComputeNodes(32),ComputeValues(32,v))
    end 
    function interpolant64(x)
        FourierInterpolantFromNodes(x,ComputeNodes(64),ComputeValues(64,v))
    end
    xs = range(0,2*pi,length=100)
    p1 = plot(xs,v.(xs),ylims=(0.3,1),color = :black,label = "")
    plot!(xs,interpolant4.(xs),linestyle = :dash,color = :blue,label ="N = 4")
    plot!(xs,interpolant8.(xs),linestyle = :dash, color = :green, label = "N = 8")
    plot!(xs,interpolant16.(xs),linestyle = :dash, color = :cyan, label = "N = 16")
    plot!(xs,interpolant32.(xs),linestyle = :dash, color = :purple, label = "N = 32")
    plot!(xs,interpolant64.(xs),linestyle = :dash, color = :red, label = "N = 64")
    #display(output)

    function interpolant4error(x)
        return abs(v(x)-FourierInterpolantFromNodes(x,ComputeNodes(4),ComputeValues(4,v)))+0.000000000000001
    end
    function interpolant8error(x)
        return abs(v(x)-FourierInterpolantFromNodes(x,ComputeNodes(8),ComputeValues(8,v)))+0.000000000000001
    end
    function interpolant16error(x)
        abs(v(x)-FourierInterpolantFromNodes(x,ComputeNodes(16),ComputeValues(16,v)))+0.000000000000001
    end
    function interpolant32error(x)
        abs(v(x)-FourierInterpolantFromNodes(x,ComputeNodes(32),ComputeValues(32,v)))+0.000000000000001
    end
    function interpolant64error(x)
        abs(v(x)-FourierInterpolantFromNodes(x,ComputeNodes(64),ComputeValues(64,v)))+0.000000000000001
    end
    xs = range(0,2*pi,length=100)
    p2 = plot(xs,interpolant4error.(xs),linestyle = :dash,color = :blue,label ="N = 4",yaxis = :log)
    plot!(xs,interpolant8error.(xs),linestyle = :dash,color = :green,label ="N = 8")
    plot!(xs,interpolant16error.(xs),linestyle = :dash,color = :cyan, label = "N = 16")
    plot!(xs,interpolant32error.(xs),linestyle = :dash,color = :purple, label = "N = 32")
    plot!(xs,interpolant64error.(xs),linestyle = :dash,color = :red, label = "N = 64")

    output = plot(p1,p2,layout =(2,1))
    #savefig(output,"FourierInterpolantPlotAndErrors.png")
end


# problem b) plot 2
let 

    function v(x)
        return 1/(2-cos(x))
    end

    function interpolant4error(x)
        return abs(v(x)-FourierInterpolantFromNodes(x,ComputeNodes(4),ComputeValues(4,v)))+0.000000000000001
    end
    function interpolant8error(x)
        return abs(v(x)-FourierInterpolantFromNodes(x,ComputeNodes(8),ComputeValues(8,v)))+0.000000000000001
    end
    function interpolant16error(x)
        abs(v(x)-FourierInterpolantFromNodes(x,ComputeNodes(16),ComputeValues(16,v)))+0.000000000000001
    end
    function interpolant32error(x)
        abs(v(x)-FourierInterpolantFromNodes(x,ComputeNodes(32),ComputeValues(32,v)))+0.000000000000001
    end
    function interpolant64error(x)
        abs(v(x)-FourierInterpolantFromNodes(x,ComputeNodes(64),ComputeValues(64,v)))+0.000000000000001
    end
    xs = range(0,2*pi,length=100)
    output = plot(xs,interpolant4error.(xs),linestyle = :dash,color = :blue,label ="N = 4",yaxis = :log)
    plot!(xs,interpolant8error.(xs),linestyle = :dash,color = :green,label ="N = 8")
    plot!(xs,interpolant16error.(xs),linestyle = :dash,color = :cyan, label = "N = 16")
    plot!(xs,interpolant32error.(xs),linestyle = :dash,color = :purple, label = "N = 32")
    plot!(xs,interpolant64error.(xs),linestyle = :dash,color = :red, label = "N = 64")
    
    display(output)
    #savefig(output,"plot4.png")
end


# problem c) plot of Lagrange polynomials.
let 
    xvalues = ComputeNodes(6)

    function lagrange0(x)
        return 1/6 * sin(3 * (x-xvalues[1]))*cot(1/2 * (x-xvalues[1]))
    end

    function lagrange1(x)
        return 1/6 * sin(3 * (x-xvalues[2]))*cot(1/2 * (x-xvalues[2]))
    end

    function lagrange2(x)
        return 1/6 * sin(3 * (x-xvalues[3]))*cot(1/2 * (x-xvalues[3]))
    end

    function lagrange3(x)
        return 1/6 * sin(3 * (x-xvalues[4]))*cot(1/2 * (x-xvalues[4]))
    end

    function lagrange4(x)
        return 1/6 * sin(3 * (x-xvalues[5]))*cot(1/2 * (x-xvalues[5]))
    end

    function lagrange5(x)
        return 1/6 * sin(3 * (x-xvalues[6]))*cot(1/2 * (x-xvalues[6]))
    end

    xs = range(0,2*pi,length = 100)

    output = plot(xs,lagrange0.(xs),label = "h₀")
    plot!(xs,lagrange1.(xs),label = "h₁")
    plot!(xs,lagrange2.(xs),label = "h₂")
    plot!(xs,lagrange3.(xs),label = "h₃")
    plot!(xs,lagrange4.(xs),label = "h₄")
    plot!(xs,lagrange5.(xs),label = "h₅")

    #savefig(output, "LagrangePolynomialPlot.png")
end




#plot of the continuous fourier coefficients. 
#=
N = 16
û = Array{Float64,1}(zeros(N))
for j in 1:(N)
    û[j] = 2/sqrt(3) * (2-sqrt(3))^j 
end

plot1 = scatter(û,xlabel = "n",label = "aₙ",title = "Continuous Fourier Coefficients",
        yticks = [10^-8,10^-6,10^-4,10^-2,1],
        yaxis = :log,
        marker = :diamond,
        color = :red)
        
=#




# problem e) plot 1
let 
    function w(x)
        return exp(sin(x))
    end

    function wderiv(x)
        return exp(sin(x))*cos(x)
    end

    wvalues4 = ComputeValues(4,w)
    wvalues8 = ComputeValues(8,w)
    wvalues16 = ComputeValues(16,w)

    derivmatrix4 = FourierDerivativeMatrix(4)
    derivmatrix8 = FourierDerivativeMatrix(8)
    derivmatrix16 = FourierDerivativeMatrix(16)
    wderivapprox4 = MxVDerivative(derivmatrix4,wvalues4,1)
    wderivapprox8 = MxVDerivative(derivmatrix8,wvalues8,1)
    wderivapprox16 = MxVDerivative(derivmatrix16,wvalues16,1)
    #uderiv2 = real(FourierDerivativeByFFT(wvalues32,1))

    xs = range(0,2*pi,length=100)

    p1 = scatter(ComputeNodes(4),wderivapprox4,marker = :cross,color = :red,label = "d/dx I₄w(x)")
    plot!(xs,wderiv.(xs),color = :black,label ="w' ")

    p2 = scatter(ComputeNodes(8),wderivapprox8,marker = :cross,color = :red,label = "d/dx I₈w(x)")
    plot!(xs,wderiv.(xs),color = :black,label ="w' ")

    p3 = scatter(ComputeNodes(16),wderivapprox16,marker = :cross,color = :red,label = "d/dx I₁₆w(x)")
    plot!(xs,wderiv.(xs),color = :black,label ="w' ")
    
    wderivvalues4 = ComputeValues(4,wderiv)
    wderivvalues8 = ComputeValues(8,wderiv)
    wderivvalues16 = ComputeValues(16,wderiv)

    p4 = plot(abs.(wderivapprox4-wderivvalues4),label = "error",color = :red)
    p5 = plot(abs.(wderivapprox8-wderivvalues8),label = "error",color = :red)
    p6 = plot(abs.(wderivapprox16-wderivvalues16),label = "error",color = :red)


    output = plot(p1,p2,p3,p4,p5,p6,layout = (2,3))
    #savefig(output,"MatrixDerivativePlots.png")
end

# code for problem f)
let 
    
    derivmatrix32 = FourierDerivativeMatrix(32)
    derivmatrix64 = FourierDerivativeMatrix(64)
    derivmatrix128 = FourierDerivativeMatrix(128)

    function w0(x)
        if( 0 <= x < pi)
            return -cos(2*x)
        elseif( pi<= x <= 2*pi)
            return cos(2*x)
        end
    end

    w0values32 = ComputeValues(32,w0)
    w0values64 = ComputeValues(64,w0)
    w0values128 = ComputeValues(128,w0)

    function w1(x)
        if( 0 <= x < pi)
            return -sin(2*x)*1/2
        elseif( pi<= x <= 2*pi)
            return sin(2*x)*1/2
        end
    end

    w1values32 = ComputeValues(32,w1)
    w1values64 = ComputeValues(64,w1)
    w1values128 = ComputeValues(128,w1)

    function w2(x)
        if( 0 <= x < pi)
            return cos(2*x)*1/4
        elseif( pi<= x <= 2*pi)
            return -cos(2*x)*1/4
        end
    end

    w2values32 = ComputeValues(32,w2)
    w2values64 = ComputeValues(64,w2)
    w2values128 = ComputeValues(128,w2)

    function w3(x)
        if( 0 <= x < pi)
            return sin(2*x)*1/8
        elseif( pi<= x <= 2*pi)
            return -sin(2*x)*1/8
         end
    end

    w3values32 = ComputeValues(32,w3)
    w3values64 = ComputeValues(64,w3)
    w3values128 = ComputeValues(128,w3)

    w0approx32 = MxVDerivative(derivmatrix32,w1values32,1)
    w0approx64 = MxVDerivative(derivmatrix64,w1values64,1)
    w0approx128 = MxVDerivative(derivmatrix128,w1values128,1)

    w1approx32 = MxVDerivative(derivmatrix32,w2values32,1)
    w1approx64 = MxVDerivative(derivmatrix64,w2values64,1)
    w1approx128 = MxVDerivative(derivmatrix128,w2values128,1)

    w2approx32 = MxVDerivative(derivmatrix32,w3values32,1)
    w2approx64 = MxVDerivative(derivmatrix64,w3values64,1)
    w2approx128 = MxVDerivative(derivmatrix128,w3values128,1)

    #=
    xs = range(0,2*pi,length=100)
    w0plot = plot(xs,w0.(xs),color = :black,label = "w⁰",thickness_scaling = 1.5)
    scatter!(ComputeNodes(32),w0approx32,marker = :cross,label = "N = 32",color = :blue)
    scatter!(ComputeNodes(64),w0approx64,marker = :cross,label = "N = 64",color = :red)
    scatter!(ComputeNodes(128),w0approx128,marker = :cross,label = "N = 128",color = :green)

    w1plot = plot(xs,w1.(xs),color = :black,label = "w¹",thickness_scaling = 1.5)
    scatter!(ComputeNodes(32),w1approx32,marker = :cross,label = "N = 32",color = :blue)
    scatter!(ComputeNodes(64),w1approx64,marker = :cross,label = "N = 64",color = :red)
    scatter!(ComputeNodes(128),w1approx128,marker = :cross,label = "N = 128",color = :green)

    w2plot = plot(xs,w2.(xs),color = :black,label = "w²",thickness_scaling = 1.5)
    scatter!(ComputeNodes(32),w2approx32,marker = :cross,label = "N = 32",color = :blue)
    scatter!(ComputeNodes(64),w2approx64,marker = :cross,label = "N = 64",color = :red)
    scatter!(ComputeNodes(128),w2approx128,marker = :cross,label = "N = 128",color = :green)

    output = plot(w0plot,w1plot,w2plot, layout = (3,1))
    =#

    w2values1024 = ComputeValues(1024,w2)
    w3values1024 = ComputeValues(1024,w3)
    derivmatrix1024 = FourierDerivativeMatrix(1024)
    w2approx1024 = MxVDerivative(derivmatrix1024,w3values1024,1)

    error = 0
    for j in 1:1024
        error = error + (w2values1024[j] - w2approx1024[j])^2
    end
    error = error * 2*pi/1024
    #savefig(output, "plot7.png")
end




# code for problem g)
let 
    
    function w(x)
         return exp(sin(x))
    end

    t = @timed MxVDerivative(FourierDerivativeMatrix(32),ComputeValues(32,w),1)
    t.time
    real(FourierDerivativeByFFT(ComputeValues(32,w),1))

    matrixTimes50 = Array{Float64,1}(zeros(50))
    for j in 1:50
        t = @timed MxVDerivative(FourierDerivativeMatrix(2*j),ComputeValues(2*j,w),1)
        matrixTimes50[j] = t.time
    end
    matrixTimes500 = Array{Float64,1}(zeros(500))
    for j in 1:500
        t = @timed MxVDerivative(FourierDerivativeMatrix(2*j),ComputeValues(2*j,w),1)
        matrixTimes500[j] = t.time
    end

    FFTTimes50 = Array{Float64,1}(zeros(50))
    for j in 1:50
          t = @timed FourierDerivativeByFFT(ComputeValues(2*j,w),1)
         FFTTimes50[j] = t.time
    end
    FFTTimes500 = Array{Float64,1}(zeros(500))
    for j in 1:500
          t = @timed FourierDerivativeByFFT(ComputeValues(2*j,w),1)
         FFTTimes500[j] = t.time
    end

    evenNumbers50 = Array{Float64,1}(zeros(50))
    for j in 1:50
        evenNumbers50[j] = 2*j
    end
    evenNumbers500 = Array{Float64,1}(zeros(500))
    for j in 1:500
        evenNumbers500[j] = 2*j
    end

    plot50 = plot(evenNumbers50,matrixTimes50,label = "Matrix based",color = :red, linestyle = :dash,ylabel = "computation time (s)",xlabel = "N")
    plot!(evenNumbers50,FFTTimes50,color = :blue, label = "FFT-based")

    plot500 = plot(evenNumbers500,matrixTimes500,label = "Matrix based",color = :red, linestyle = :dash,ylabel = "computation time (s)",xlabel = "N")
    plot!(evenNumbers500,FFTTimes500,color = :blue, label = "FFT-based")

    output = plot(plot50,plot500,layout = (2,1))
    #savefig(output,"plot8.png")
end


# code for problem h)
let 
    xs = range(-1,1,length=100)
    p1 = plot(xs,LegendrePolynomial.(0,xs),label = "L₀")
    plot!(xs,LegendrePolynomial.(1,xs),label = "L₁")
    plot!(xs,LegendrePolynomial.(2,xs),label = "L₂")
    plot!(xs,LegendrePolynomial.(3,xs),label = "L₃")
    plot!(xs,LegendrePolynomial.(4,xs),label = "L₄")
    plot!(xs,LegendrePolynomial.(5,xs),label = "L₅")

    p2 = plot(xs,ChebyshevPolynomial.(0,xs),label = "T₀")
    plot!(xs,ChebyshevPolynomial.(1,xs),label = "T₁")
    plot!(xs,ChebyshevPolynomial.(2,xs),label = "T₂")
    plot!(xs,ChebyshevPolynomial.(3,xs),label = "T₃")
    plot!(xs,ChebyshevPolynomial.(4,xs),label = "T₄")
    plot!(xs,ChebyshevPolynomial.(5,xs),label = "T₅")

    output = plot(p1,p2,layout = (2,1))
    #savefig(output,"LegendreChebyshevPlots.png")
end

# More code for problem i) Normalized Legendre polynomials.

let 
    xs = range(-1,1,length=100)
    output = plot(xs,OrthonormalJacobiPolynomial.(0.0,0.0,0,xs),label = "̃L₀")
    plot!(xs,OrthonormalJacobiPolynomial.(0.0,0.0,1,xs),label = "̃L₁")
    plot!(xs,OrthonormalJacobiPolynomial.(0.0,0.0,2,xs),label = "̃L₂")
    plot!(xs,OrthonormalJacobiPolynomial.(0.0,0.0,3,xs),label = "̃L₃")
    plot!(xs,OrthonormalJacobiPolynomial.(0.0,0.0,4,xs),label = "̃L₄")
    plot!(xs,OrthonormalJacobiPolynomial.(0.0,0.0,5,xs),label = "̃L₅")

    #savefig(output,"NormalizedLegendrePlot.png")
end




# code for problem i)
let 

    function u(x)
        return 1/(2-cos(2*pi*x))
    end
    #xs = range(-1,0.9,length=100)

    xvalues10,wvalues10 = gausslegendre(10)
    xvalues40,wvalues40 = gausslegendre(40)
    xvalues80,wvalues80 = gausslegendre(80) 
    xvalues100,wvalues100 = gausslegendre(100)
    xvalues200,wvalues200 = gausslegendre(200)

    coeffs10 = DiscreteLegendreCoefficients(200,xvalues10,wvalues10,u)
    coeffs40 = DiscreteLegendreCoefficients(200,xvalues40,wvalues40,u)
    coeffs80 = DiscreteLegendreCoefficients(200,xvalues80,wvalues80,u)
    coeffs100 = DiscreteLegendreCoefficients(200,xvalues100,wvalues100,u)
    coeffs200 = DiscreteLegendreCoefficients(200,xvalues200,wvalues200,u)


    p1 = scatter(abs.(coeffs10),ylims =(10^-19,1),yaxis = :log, marker = :cross, label = "N = 10")
    scatter!(abs.(coeffs40),ylims =(10^-19,1),yaxis = :log, marker = :cross,label = "N = 40")
    scatter!(abs.(coeffs80),ylims =(10^-19,1),yaxis = :log, marker = :cross, label = "N = 80")


    p2 = scatter(abs.(coeffs100),ylims =(10^-19,1),yaxis = :log, marker = :cross,label = "N = 100")
    scatter!(abs.(coeffs200),ylims =(10^-19,1),yaxis = :log, marker = :cross,label = "N = 200")

    output = plot(p1,p2,layout = (2,1))
    #savefig(output, "LegendreCoefficients.png")
end


# code for problem j)
let 
    xv,w = gausslobatto(6)

    function interpolant0(x)
        return LagrangeInterpolant(0.0,0.0,xv,x)[1]
    end

    function interpolant1(x)
         return LagrangeInterpolant(0.0,0.0,xv,x)[2]
    end

    function interpolant2(x)
          return LagrangeInterpolant(0.0,0.0,xv,x)[3]
    end

    function interpolant3(x)
           return LagrangeInterpolant(0.0,0.0,xv,x)[4]
    end

    function interpolant4(x)
           return LagrangeInterpolant(0.0,0.0,xv,x)[5]
    end

    function interpolant5(x)
        return LagrangeInterpolant(0.0,0.0,xv,x)[6]
    end

    xs = range(-1,1,length = 100)


    output = plot(xs,interpolant0.(xs),label = "h₀")
    plot!(xs,interpolant1.(xs),label = "h₁")
    plot!(xs,interpolant2.(xs),label = "h₂")
    plot!(xs,interpolant3.(xs),label = "h₃")
    plot!(xs,interpolant4.(xs),label = "h₄")
    plot!(xs,interpolant5.(xs),label = "h₅")

    #savefig(output,"LagrangeInterpolants.png")
end


# code for problem k) Plots are currently marked as comments.

let 
    function v(x)
        return exp(-sin(pi*x))
    end
    
    x10,w10 = gausslobatto(10)
    V10 = VandermondeMatrix(0.0,0.0,x10)
    Vx10 = VandermondeDerivativeMatrix(0.0,0.0,x10)
    D10 = Vx10 * inv(V10)
    vapprox10 = D10 * v.(x10)

    x40,w40 = gausslobatto(40)
    V40 = VandermondeMatrix(0.0,0.0,x40)
    Vx40 = VandermondeDerivativeMatrix(0.0,0.0,x40)
    D40 = Vx40 * inv(V40)
    return vapprox40 = D40 * v.(x40)
    
    x80,w80 = gausslobatto(80)
    V80 = VandermondeMatrix(0.0,0.0,x80)
    Vx80 = VandermondeDerivativeMatrix(0.0,0.0,x80)
    D80 = Vx80 * inv(V80)
    vapprox80 = D80 * v.(x80)
    
    function vderiv(x)
        return -pi * exp(-sin(pi*x))*cos(pi*x)
    end
    
    
    xs = range(-1,1,length = 100)

    p1 = plot(xs,vderiv.(xs),label = "v'")
    scatter!(x10,vapprox10,marker = :cross,label = "N = 10")

    p4 = plot(abs.(vderiv.(x10) - vapprox10),label = "error",color = :red)

    p2 = plot(xs,vderiv.(xs),label = "v'")
    scatter!(x40,vapprox40,marker = :cross, label = "N = 40")

    p5 = plot(abs.(vderiv.(x40) - vapprox40),label = "error",color = :red)

    p3 = plot(xs,vderiv.(xs),label = "v'")
    scatter!(x80,vapprox80,marker = :cross, label = "N = 80")

    p6 = plot(abs.(vderiv.(x80) - vapprox80),label = "error",color = :red)

    output = plot(p1,p2,p3,p4,p5,p6, layout = (2,3))
    
    #savefig(output,"LegendreDerivativeApproximation.png")

    
    M10 = inv((V10) * transpose(V10))
    norm10 = sum(vapprox10.*(M10 * vapprox10))
    #normt10 = sum(vderiv.(x10).*(M10 * vderiv.(x10)))

    #M40 = inv((V40)*transpose(V40))
    #norm40 = sum(vapprox40.*(M40 * vapprox40))
    #normt40 = sum(vderiv.(x40).*(M40 * vderiv.(x40)))

    #M80 = inv((V80)*transpose(V80))
    #norm80 = sum(vapprox80.*(M80 * vapprox80))
    #normt80 = sum(vderiv.(x80).*(M80 * vderiv.(x80)))
    
end

# problem label) Not commented.
let 
    x40,w40 = gausslobatto(40)
    x80,w80 = gausslobatto(80)
    V40 = VandermondeMatrix(0.0,0.0,x40)
    V80 = VandermondeMatrix(0.0,0.0,x80)
    M40 = inv((V40)*transpose(V40))
    M80 = inv((V80)*transpose(V80))
    
    ones = Array{Float64,1}(fill(1.0,40))
    
    sum(ones.*(M40*ones))
    
    function gsin(x)
        return sin(x-1)
    end
    
    gsin40 = gsin.(x40)
    gsin80 = gsin.(x80)
    
    
    sum(gsin40.*(M40*gsin40))
    sum(gsin80.*(M80*gsin80))
    
    V40*transpose(V40)
    
end 












 