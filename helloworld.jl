using Plots

# println("Hello :)")


function my_func(x)
    return x.*x
end

function make_plot()
    # x = LinRange(0,10, 100)
    # y = my_func(x)

    x = [1, 2, 3]
    y = [1, 2, 3]

    # g = plot(x,y)
    # display(g)
    p = plot(x,y)
    # xlabel("X")
    # ylabel("Y")
    # PyPlot.title("Your Title Goes Here")
    # grid("on")
end

if abspath(PROGRAM_FILE) == @__FILE__
    # A = range(start = 1, length=10, stop=100)
    # B = LinRange(1, 5.5, 9)
   
    
    # make_plot()
    
    # println(length(B))
    # for i in B
    #     println(i)
    # end
    # println(B)

    data = rand(10)
    Plots.pyplot()
    f = plot(data);
    display(f)
    readline()
    # @show f
    # g = plot(x,y)
    # display(g)
    # println(typeof(A))
    # y = multiply(A, B)
    # println(y)
end








