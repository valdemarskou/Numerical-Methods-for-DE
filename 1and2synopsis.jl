
#using CairoMakie
using Plots


# function is abs(x). Truncation is P_2N f
function absError(N)
    τ = (pi^2)/12
    for k in 1:Int(2*N)
        τ = τ - 4*(1-(-1)^k)/(pi^2 * k^4)
    end

    return τ
end

# function is 3/(5-4cos(x)). Truncation is again P_2N f
function f2Error(N)
    return 2/(3* 4^(2*N))    
end

#function is sum_n=1^inf sin(nx)/n. Truncation is again P_2N f1
function f3Error(N)
    τ = (pi^2)/12
    for k in 1:Int(2*N)
        τ = τ - 1/(2*k^2)
    end

    return τ
end

let 
    xs = range(1,length =16)
    ys1 = absError.(xs)
    ys2 = f2Error.(xs)
    ys3 = f3Error.(xs)

    plot = scatter(xs, ys1,
        yaxis = :log,
        #yticks = [10^-5,10^-4,10^-3,10^-2],
        marker = :diamond,
        color = :red,
        label = "f₁",
        title = "Truncation errors L^2 norm")

    scatter!(xs, ys2,
        yaxis = :log,
        #yticks = [10^-5,10^-4,10^-3,10^-2],
        marker = :circle,
        color = :blue,
        label = "f₂",
        title = "Truncation errors L^2 norm",
        )

    scatter!(xs, ys3,
        yaxis = :log,
        #yticks = [10^-5,10^-4,10^-3,10^-2],
        marker = :rect,
        color = :yellow,
        label = "f₃",
        title = "Truncation errors L^2 norm",
        )
#savefig(plot,"truncationerrors.png")
display(plot)    
end



#blah blah blah
function f₂(x)
    return 3/(5-4*cos(x))
end


using Plots

plot2 = plot(f₂, 0:0.1:2*pi,color = :red, label = "f₂(x)",
    xticks = [0,pi/2,pi,3*pi/2,2*pi],
    formatter = :plain,
    xlims=(0,2*pi))

#Adding the N projections to the plot.


    

for N in [2,4,8,16,32,64]
    function f2trunc(x)
        val = 0
        for k in 1:N
            val = val + 1/(2^(k)) * cos(k*x)
        end
        return val
    end

    plot!(f2trunc, 0:0.1:2*pi,color = :red, label = "f₂(x)",
    xticks = [0,pi/2,pi,3*pi/2,2*pi],
    formatter = :dot,
    xlims=(0,2*pi))



end


# plot!()