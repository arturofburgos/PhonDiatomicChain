#========================================#
# Name: Arturo Burgos                    #
#                                        #
#                                        #
#           Spring Mass System           #
#                                        #
#========================================#

# Second order ODE: Harmonic Oscillation - one spring and one mass

# Create the time grid: 

T = 20 # Lenght of the time interval
dt = 0.002 # Time-step size
N = Int(T/dt) + 1 # Number of time-steps
t = range(0, T, N+1)

# Constants: 

k = 1.9802
m = 0.0689
c = 0.05


M = [1 0; 0 m]
K = [0 -1; k 0]
F = [0; 0.1]

C = [0 0; 0 c]

I = [1 0; 0 1]




# Using Eulers Method: 

# Set initial conditions:

x0 = 0.0 # Initial Position 
ẋ0 = 0.0 # Initial Velocity  

u0 = [x0, ẋ0] # Initial condition matrix 



# Array to store the position at each time-step: 


function trapezoidal(y0, h, N, M, C, K, F, I)

    save_var = [copy(y0[1])]
    y = copy(y0)

    for i in 1:N
        y[:] = inv(I + (h/2) * (inv(M)*C + inv(M)*K)) * (I - (h/2) * (inv(M)*C + inv(M)*K)) * y[:] +
        inv(I + (h/2) * (inv(M)*C + inv(M)*K)) * inv(M) * F * h # Or h only it is h only
        push!(save_var,copy(y[1]))
    end
    return save_var
end

x = trapezoidal(u0, dt, N, M, C, K, F, I)


using Plots, LaTeXStrings, Measures # Possible to use PyPlot too.
gr()


@userplot springmass
@recipe function f!(var::springmass)
        
    size --> (950,400)
    margin --> 5mm
    top_margin --> 3mm
    framestyle --> :box
    grid --> :true
    gridalpha --> 5
    gridwidth --> 0.3
    minorgrid --> :true
    minorgridalpha --> 1
    minorgridwidth --> 0.05
    fontfamily --> "Computer Modern"

    t, x = var.args
    title --> "Mass position over time"
    legend --> :topright
    
    
    xaxis --> ("Time [s]")
    yaxis -->("Position [m]")
    #yguidefontrotation --> -90
    seriestype --> :line 
    linewidth --> 1.5
    #color --> :red 
    
    return t, x

end


figure = springmass(t,x[:,1], label = "Numerical")
display(figure)

max_disp = maximum(x[:,1])

# Compare with the Analytical Solution: 

#x_exact = x0 * cos.(c*t) + (ẋ0/c) * sin.(c*t)
#springmass!(t,x_exact, label = "Analytical")