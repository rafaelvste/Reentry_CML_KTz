# Eckmann-Ruelle KTz model LxL network

#Simulates the lattice of diffusively coupled KTz cells and calculates the Lyapunov spectrum.
#The exponents are approximated using the Eckmann-Ruelle method, where the Jacobian
#matrix is triangularized by LU decompostion at each iteration using SciPy LU function.
#Outputs the main exponent against coupling J.

# Import the libraries:
import numpy as np
import random as rd
import matplotlib.pyplot as plt
from scipy.linalg import lu
import time

# Starting time
start = time.time()

# Random seed:
rd.seed(666)

# Parameters for the KTz model:			  			  
k = 0.6
T = 0.154
delt = 0.001
lamb = 0.001
Xr = -0.48
Ie = 0.0

# Diffusive coupling constant:		  
J_0 = 0.01
J_f = 0.04
dJ = 1.0e-4

# Network parameters	  
L = 2
N = L**2
dt = 92 #Pacing period (P in the paper)

# Discretize J variable and store values in a tensor:
n_J_max = int((J_f - J_0) / dJ) + 1
J = J_0 - dJ
J_val = np.zeros((N,4,n_J_max))
values_J = np.zeros((n_J_max,1))
for n_J in range(0, n_J_max): #1:n_J_max
    J += dJ;
    J_val[:,:,n_J] = J
    # Output:
    values_J[n_J,0] = J
    
#Total time steps and transient:
t_max = 100000
t_trans = 1

# Allocate memory:
rand_y = np.zeros((N,n_J_max))
rand_x = np.zeros((N,n_J_max))
Xo = np.zeros((N,n_J_max))
Yo = np.zeros((N,n_J_max))
Zo = np.zeros((N,n_J_max))      
X = np.zeros((N,n_J_max))
Y = np.zeros((N,n_J_max))
Z = np.zeros((N,n_J_max))
I = np.zeros((N,n_J_max))
I_stim = np.zeros((N,n_J_max))
arg = np.zeros((N,n_J_max))
O = np.zeros((N,4,n_J_max))
eps = np.zeros((N,4,n_J_max))
G = np.zeros((N,4,n_J_max))
C = np.zeros((N,4,n_J_max))
# Lyapunov variables:
a0 = np.zeros((N,n_J_max))
a11 = np.zeros((N,n_J_max))
a12 = np.zeros((N,n_J_max))
a13 = np.zeros((N,n_J_max))
b11 = np.zeros((N,n_J_max))
M_ii = np.zeros((3,3,N,n_J_max))
M_ij = np.zeros((3,3,N,n_J_max))
M_0 = np.zeros((3,3,n_J_max))
M = np.zeros((3*N,3*N,n_J_max))
M_lu = np.zeros((3*N,3*N))
O1 = np.zeros((3*N,3*N,n_J_max))
O1_lu = np.zeros((3*N,3*N))
O2_lu = np.zeros((3*N,3*N))
lambda0 = np.zeros((3*N,n_J_max)) 
# Output:
exp_lyap = np.zeros((n_J_max,3*N))

# Generate random numbers for the IC and the stimulus vector (cells of the rightmost column):
for a in range(1, L+1): #1:L
    rand_y[a*L-1,:] = (rd.random() - 1.0/2.0)*1.0e-7
    rand_x[a*L-1,:] = (rd.random() - 1.0/2.0)*1.0e-7
    I_stim[a*L-1,:] = 0.1
    # Use to stimulate more columns (count from right to left):
#    for a_o in range(1, 5): #1:4
#       I[a*L-a_o,0] = 0.1
#       rand_y[a*L-1-a_o,1] = (rd.random() - 1.0/2.0)*1.0e-7
#       rand_x[a*L-1-a_o,1] = (rd.random() - 1.0/2.0)*1.0e-7     

# Number of synapse for each cell:
n_s = np.zeros((N,1))
if L == 2:
    n_s[:,0] = 2.0
if L > 2:
    n_s[0,0] = 2.0                  #1
    for a in range(1, L-1):         #2:(L-1)
        n_s[a,0] = 3.0
    n_s[L-1,0] = 2.0                #L
    n_s[L,0] = 3.0                  #L+1
    for a in range(L+1, N-L-1):     #(L+2):(N-L-1)
        n_s[a,0] = 4.0
    n_s[N-L-1,0] = 3.0              #N-L
    n_s[N-L,0] = 2.0                #N-L+1    
    for a in range(N-L+1, N-1): #(N-L+2):(N-1)
        n_s[a,0] = 3.0
    n_s[N-1,0] = 2.0                #N
# Fix the number of synapses of the left and right columns:
if L > 3:    
    for a in range(2, L-1):         #2:(L-2)
        n_s[a*L-1,0] = 3.0          #a*L
        n_s[a*L,0] = 3.0            #a*L+1

# Iterate time:
for t1 in range(1, t_max+1):    #1:t_max
    # Initial conditions for the KTz model:			 
    if t1 == 1:
        St = 1
        Zo[:,:] = 0.0
        Yo[:,:] = -0.5
        Xo[:,:] = -0.5
        # Random IC and stimulus:
        Yo += rand_y
        Xo += rand_x
        I = I_stim        					 
    # Update variables of the KTz model:
    elif t1 > 1:
        St = St + 1
        Zo = Z
        Yo = Y
        Xo = X
        # Periodic stimulus:
        if St <= 10: #Adjust the duration of the stimulus
            I += I_stim
        if St == dt:
            St = 0
    # Calculate the KTz equations:			 
    arg = (Xo - k * Yo + Zo + Ie + I) / T
    X = np.divide(arg, (1.0 + abs(arg))) #Logistic KTz
    Y = Xo
    Z = (1.0 - delt) * Zo - lamb * (Xo - Xr)
    # Lyapunov exponent (vectorized):
    if t1 >= t_trans:       
        a0 = np.divide(1.0, T*(1.0 + abs(arg))**2)
        # Values in the single element Jacobian matrices M_ii (vectorized):
        a11 = (1.0 - n_s*J_val[:,0,:]) * a0
        a12 = -k * a0
        a13 = a0
        # Single value in the cross elements Jacobian matrices M_ij (vectorized):
        b11 = J_val[:,0,:] * a0
        # Place the values in the M_ii and M_ij tensors (the 3rd dimension refer to network cell, 4th refer to J value):
        M_ii[0,0,:,:] = a11[:,:]        #1,1,:
        M_ii[0,1,:,:] = a12[:,:]        #1,2,:
        M_ii[0,2,:,:] = a13[:,:]        #1,3,:
        M_ii[1,0,:,:] = 1.0             #2,1,:
        M_ii[2,0,:,:] = -lamb           #3,1,:
        M_ii[2,2,:,:] = 1.0 - delt      #3,3,:
        M_ij[0,0,:,:] = b11[:,:]        #1,1,:
        # Fill the Jacobian matrix of the network with M_ii and M_ij:
        u_o = -3
        for u in range(0, N): #1:N
            u_o += 3
            v_o = -3
            for v in range(0, N): #1:N
                v_o += 3
                if v == u:
                    M[u_o:u_o+3, v_o:v_o+3, :] = M_ii[:,:,u,:]
                elif v == (u-1) or v == (u+1):
                    M[u_o:u_o+3, v_o:v_o+3, :] = M_ij[:,:,u,:]
                elif v == (u-L) or v == (u+L):
                    M[u_o:u_o+3, v_o:v_o+3, :] = M_ij[:,:,u,:]
                else:
                    M[u_o:u_o+3, v_o:v_o+3, :] = M_0[:,:,:]
        # Fix the matrices of the left and right columns:
        u_f = 0
        for u in range(1, L): #1:(L-1)
            u_o = u_f + 3*(L-1)
            u_f = u_o + 3
            M[u_o:u_f, u_f:u_f+3, :] = M_0[:,:,:]
            M[u_f:u_f+3, u_o:u_f, :] = M_0[:,:,:]
        # Perform the LU decomposition:
        for n_J in range(0, n_J_max): #1:n_J_max
            M_lu[:,:] = M[:,:,n_J]
            if t1 > t_trans:
                O1_lu[:,:] = O1[:,:,n_J]
                M_lu = np.matmul(M_lu, O1_lu)
    #        Per, LT, O2 = lu(M)
    #        O1 = np.matmul(Per, LT)
            O1_lu, O2_lu = lu(M_lu, permute_l=True)
            O1[:,:,n_J] = O1_lu[:,:]
            lambda_n_J = np.log(abs(np.diag(O2_lu).reshape((3*N,1))))
            lambda0[:,n_J] += lambda_n_J[:,0]       
    # Matrix filled with the differences between potentials of first nearest neighbors:
    O[0,0,:] = X[1,:] - X[0,:]            #O(1,1) = X(2,1) - X(1,1)
    O[0,1,:] = 0.0                        #O(1,2) = 0.0
    O[0,2,:] = X[L,:] - X[0,:]            #O(1,3) = X(1+L,1) - X(1,1)
    O[0,3,:] = 0.0                        #O(1,4) = 0.0
    for a in range(1, L):               #2:L
        O[a,0,:] = X[a+1,:] - X[a,:]      #O(a,1) = X(a+1,1) - X(a,1)
        O[a,1,:] = X[a-1,:] - X[a,:]      #O(a,2) = X(a-1,1) - X(a,1)
        O[a,2,:] = X[a+L,:] - X[a,:]      #O(a,3) = X(a+L,1) - X(a,1)
        O[a,3,:] = 0.0                    #O(a,4) = 0.0
    for a in range(L, N-L):             #(L+1):(N-L)
        O[a,0,:] = X[a+1,:] - X[a,:]      #O(a,1) = X(a+1,1) - X(a,1)
        O[a,1,:] = X[a-1,:] - X[a,:]      #O(a,2) = X(a-1,1) - X(a,1)
        O[a,2,:] = X[a+L,:] - X[a,:]      #O(a,3) = X(a+L,1) - X(a,1)
        O[a,3,:] = X[a-L,:] - X[a,:]      #O(a,4) = X(a-L,1) - X(a,1)
    for a in range(N-L, N-1):           #((N-L)+1):(N-1)
        O[a,0,:] = X[a+1,:] - X[a,:]      #O(a,1) = X(a+1,1) - X(a,1)
        O[a,1,:] = X[a-1,:] - X[a,:]      #O(a,2) = X(a-1,1) - X(a,1)
        O[a,2,:] = 0.0                    #O(a,3) = 0.0
        O[a,3,:] = X[a-L,:] - X[a,:]      #O(a,4) = X(N-L,1) - X(N,1)    
    O[N-1,0,:] = 0.0                      #O(N,1) = 0.0
    O[N-1,1,:] = X[N-2,:] - X[N-1,:]      #O(N,2) = X(N-1,1) - X(N,1)
    O[N-1,2,:] = 0.0                      #O(N,3) = 0.0
    O[N-1,3,:] = X[N-1-L,:] - X[N-1,:]    #O(N,4) = X(N-L,1) - X(N,1)             
    # Fix the boundary conditions of the left and right columns:
    for a in range(1, L):               #1:(L-1)
        O[a*L-1,0,:] = 0.0                #O(a*L,1) = 0.0
        O[a*L,1,:] = 0.0                  #O(a*L+1,2) = 0.0
    # Diffusive coupling:
#    eps[:,:] = 0.0
    G = J_val #* (1.0 + eps*1.0e-5)
    C = G * O
    I = np.sum(C,axis=1)
    
exp_lyap = lambda0.T / (t_max - t_trans)

# End time
end = time.time()

# Plotting the results:
plt.plot(values_J, exp_lyap[:,0], '-o', markersize = 2, linewidth = 1, color = 'blue')
plt.xlabel('J')
plt.ylabel('Lyapunov exponent')
plt.grid()
plt.show()

# Save to file
data = np.concatenate((values_J, exp_lyap), axis=1)
np.savetxt("exp_lyap_python.dat", data, delimiter=" ")

data = np.array((values_J[:,0], exp_lyap[:,0])).T
np.savetxt("exp_lyap1_python.dat", data, delimiter=" ")

run_time = end-start
runtime = open('runtime.txt', 'w')
runtime.write(str(run_time))
runtime.close()