

using FFTW
include("DiscreteFourierTransform.jl")

# Testing documentation of the FFTW package
#=
N=2^9
test = Array{ComplexF64,1}([(exp(-2*pi*im/N))^2,exp(-2*2*pi*im/N)])

testInput = Array{ComplexF64,1}(sin.(1:2^9)+cos.(1:2^9)+(1:2^9)/(2^9))
FFTWOutput = fft(testInput)


# Plot of the coefficients
scatter(FFTWOutput,title = "Fourier Transform",xlabel = "Real Part", ylabel = "Imaginary Part", label = false, color =:blue)
=#


# FFFT of two real vectors. Seems to work fine but I am confused about the indexing.

function FFFTOfTwoRealVectors(xvalues::Array{Float64,1},yvalues::Array{Float64,1})
    @assert length(xvalues) == length(yvalues)
    N = length(xvalues)
    zvalues = Array{ComplexF64,1}()
    for k=1:N
        push!(zvalues,xvalues[k]+im*yvalues[k])
    end
    ztransform = fft(zvalues)

    xtransform = Array{ComplexF64,1}()
    push!(xtransform,real(ztransform[1]))
    ytransform = Array{ComplexF64,1}()
    push!(ytransform,imag(ztransform[1]))

    for k=2:N 
        push!(xtransform,(ztransform[k]+conj(ztransform[N+2-k]))/(2))
        push!(ytransform, -im*(ztransform[k]-conj(ztransform[N+2-k]))/(2))
    end

    return xtransform, ytransform
end


# BFFT for two real vectors. Total hack-job because of indexing, using both push and fill.

function BFFTForTwoRealVectors(Xvalues::Array{Float64,1},Yvalues::Array{Float64,1})
    @assert length(Xvalues) == length(Yvalues)
    N = length(Xvalues)
    Zvalues = Array{ComplexF64,1}()
    for k=1:N 
        push!(Zvalues, Xvalues[k] + im*Yvalues[k])
    end
    Zvalues = bfft(Zvalues)

    xvalues = Array{ComplexF64,1}(zeros(N))
    xvalues[1] = real(Zvalues[1])
    yvalues = Array{ComplexF64,1}(zeros(N))
    yvalues[1] = imag(Zvalues[1])

    for k=1:(N-1)
        xvalues[k+1] = (Zvalues[k+1]+conj(Zvalues[N+1-k]))/(2)
        yvalues[k+1] = -im*(Zvalues[k+1]-conj(Zvalues[N+1-k]))/(2)
    end

    return xvalues, yvalues
end


#Foward DFT by even-odd decomposition. Requires N to be even. There is seemingly an
# error with the real parts of the output, unsure how to fix.

# # # HAS ERRORS # # # 

function FFFTEO(fvalues::Array{Float64,1},wvalues::Array{ComplexF64,1})
    N = length(fvalues)
    Zvalues = Array{ComplexF64,1}(zeros(Int(N/2)))
    
    for j=0:(Int(N/2) - 1) 
        Zvalues[j+1] = fvalues[2*j+1] + im*fvalues[2*j+2]
        # Dette er rigtigt.
    end

    Zvalues = fft(Zvalues)

    Fvalues = Array{ComplexF64,1}(zeros(N))
    Fvalues[1] = (real(Zvalues[1]) + imag(Zvalues[1]))
    Fvalues[Int(N/2)+1] = (real(Zvalues[1]) - imag(Zvalues[1]))

    for k=1:(Int(N/2) - 1)
        Fvalues[k+1] = ((Zvalues[k+1]+conj(Zvalues[Int(N/2)-k+1]))-im*wvalues[k+1]*(Zvalues[k+1]-conj(Zvalues[Int(N/2)-k+1])))/(2)
    end
    for k=1:(Int(N/2) - 1) 
        Fvalues[N-k+1] = conj(Fvalues[k+1])
    end

    return Fvalues
end


# The Forward Real FFT (does not use the FFFTEO algorithm because of errors).
# Algorithm 15. Requires N to be even.

function ForwardRealFFT(xvalues::Array{Float64,1})
    N = length(xvalues)
    Xvalues = fft(xvalues)

    avalues = Array{Float64,1}(zeros(Int(N/2)+1))
    bvalues = Array{Float64,1}(zeros(Int(N/2)+1))

    for k in 0:Int(N/2)
        avalues[k+1] = 2*real(Xvalues[k+1])
        bvalues[k+1] = -2*imag(Xvalues[k+1])
    end
    bvalues[1] = 0
    bvalues[Int(N/2)+1] = 0
    
    return avalues, bvalues
end


#
# Algorithm 17

function FourierDerivativeByFFT(f::Array{Float64,1},m::Int)
   
    N = length(f) 
    F = fft(f)

    # Requires N to be even.

    for k in 0:(Int(N/2) - 1)
        F[k+1] = (im*k)^(m) * F[k+1]
    end

    if iseven(m) == false
        F[Int(N/2) + 1] = 0
    end

    for k in (Int(N/2) + 1):(N-1)
        F[k+1] = (im*(k-N))^m * F[k+1]
    end

    Df = bfft(F)/N

    
end

# Computation of the Fourier Derivative Matrix Using the Negative Sum Trick
# Algorithm 18

function FourierDerivativeMatrix(N::Int)
    
    D = Array{Float64,2}(fill(0,(N,N)))
    for i in 0:(N-1) 
        for j in 0:(N-1)
            if j != i
                D[i+1,j+1] = 1/2 * (-1)^(i-j) * cot(pi*(i-j)/N)
                D[i+1,i+1] = D[i+1,i+1] - D[i+1,j+1]
            end
        end
    end

    return D
end


# A Matrix-Vector Multiplication Procedure.
# Algorithm 19.

function MxVDerivative(D::Array{Float64,2},f::Array{Float64,1},s::Int)
    e = length(f)
    @assert((e,e) == size(D))

    InterpolantDerivative = Array{Float64,1}(zeros(e-s+1))

    for i in s:e 
        t = 0
        for j in s:e
            t = t + D[i,j]*f[j]
        end
        InterpolantDerivative[i-s+1] = t
    end

    return InterpolantDerivative    
end












# Simple testing
#=
realVector1 = Array{Float64,1}(sin.(1:2^9)+cos.(1:2^9)+(1:2^9)/(2^9))
realVector2 = Array{Float64,1}(cos.(1:2^9)+cos.(1:2^9)+sin.(1:2^9))

BFFTForTwoRealVectors(realVector1, realVector2)

wvalues = InitializeFFT(2^9,1)
FFTEOOutput = FFFTEO(realVector1,wvalues)

testOutput = fft(realVector1)

scatter(FFTEOOutput)
scatter(testOutput)
=#