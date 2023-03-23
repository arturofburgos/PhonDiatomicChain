#========================================#
# Name: Arturo Burgos                    #
#                                        #
#           Spring Mass System           #
#                                        #
#========================================#

# Second order ODE: Harmonic Oscillation - nm spring and nm mass
using LinearAlgebra, LinearMaps, Plots
#============#
# USER INPUT #
#============#

# Here define the number of pair spring-mass: 
ncells = 5 

nm = 2*ncells


# Assuming all the masses are equal:
m1 = 0.01

m2 = 0.02

m = [m1, m2]
# Assuming all springs are equal:
k1 = 22

k2 = 18 

k = [k1, k2]


# Define A and At matrices: We either can create the actual A matrix
# Or we can apply a function Atimes! that maps a matrix. 

#========#
# Matrix #
#========#

A = zeros(nm,nm)

for i in 1:nm
    for j in 1:nm
        if j == i
            A[i,j] = 1
        end
        if j == i - 1
            A[i,j] = -1
        end
    end 
end

At = A'


# Define K̃ matrix
k_vec = repeat(k, ncells)
K̃ = diagm(k_vec)

# Define K̂
K = A'*K̃*A

F = zeros(nm)
F[end] = 0.1
 
x = inv(K)*F