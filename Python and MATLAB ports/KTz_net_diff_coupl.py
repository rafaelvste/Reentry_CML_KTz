# -*- coding: utf-8 -*-

# Simulates the lattice of diffusively coupled KTz cells and returns the
# network actvity (number of cells with x > 0), membrane potential and
# current of one of the cells.

# Import the libraries:
import random as rd
import numpy as np
import matplotlib.pyplot as plt

#KTz model
def ktz_model(N, X_o, Y_o, Z_o, I_t):
    # Create array:
    arg = np.zeros((N,1))
    # Model paramaters
    k = 0.6
    T = 0.154
    delt = 0.001
    lamb = 0.001
    Xr = -0.48
    Ie = 0.0
    # Model equations:
    arg = (X_o - k * Y_o + Z_o + Ie + I_t) / T
    X_f = np.divide(arg, (1.0 + abs(arg))) #Logistic KTz
#    X_f = np.tanh(arg) #Hyperbolic KTz
    Y_f = X_o
    Z_f = (1.0 - delt) * Z_o - lamb * (X_o - Xr)
    return(X_f, Y_f, Z_f)

#Diffusive coupling (electrical synapse)
def diff_coupl(L, N, X_s):
    # Create arrays:
    O = np.zeros((N,4))
    G = np.zeros((N,4))
    C = np.zeros((N,4))
    # Diffusive coupling constant:
    J = 0.0257
    # Matrix filled with the differences between potentials of first nearest neighbors
    for a in range(0, N): #Numbering from left to right and from bottom to top
        # Connection with the right neighbor:
        if a < (N-1):                           #a < N
            O[a,0] = X_s[a+1,0] - X_s[a,0]
        # Connection with the left neighbor:
        if a > 0:                               #a > 1
            O[a,1] = X_s[a-1,0] - X_s[a,0]
        # Connection with the neighbor above:
        if a <= (N-L-1):                        #a <= (N-L)
            O[a,2] = X_s[a+L,0] - X_s[a,0]
        # Connection with the neighbor below:
        if a > (L-1):                           #a > L
            O[a,3] = X_s[a-L,0] - X_s[a,0]
    # Fix the remaining boundary conditions of the left and right columns:
    for a in range(1, L):                       #1:(L-1)
        O[a*L-1,0] = 0.0                        #O(a*L,1) = 0.0
        O[a*L,1] = 0.0                          #O(a*L+1,2) = 0.0
    # Diffusive coupling equation:
    G[:,:] = J
    C = G * O
    I_s = np.sum(C,axis=1).reshape((N,1))
    return I_s

# Random seed:
rd.seed(666)

# Network parameters
L_side = 80
N_cells = L_side**2

# Cell for measures:
cell_sample = int((N_cells-L_side)/2)

# Additional parameters:
t_max = 50000 #Total time steps
dt = 92 #Pacing period (P in the paper)
stim_dur = 10 #Stimulus duration

# Create arrays:
rand_y = np.zeros((N_cells,1))
rand_x = np.zeros((N_cells,1))
Xo = np.zeros((N_cells,1))
Yo = np.zeros((N_cells,1))
Zo = np.zeros((N_cells,1))
X = np.zeros((N_cells,1))
Y = np.zeros((N_cells,1))
Z = np.zeros((N_cells,1))
I = np.zeros((N_cells,1))
I_stim = np.zeros((N_cells,1))
I_snps = np.zeros((N_cells,1))

# Outputs:
s = np.zeros((t_max,1))
time = np.zeros((t_max,1))
potential = np.zeros((t_max,1))
current = np.zeros((t_max,1))

# Generate random numbers for the IC and the stimulus vector (cells of the rightmost column):
for u in range(1, L_side+1): #1:L_side
    rand_y[u*L_side-1,0] = (rd.random() - 1.0/2.0)*1.0e-7
    rand_x[u*L_side-1,0] = (rd.random() - 1.0/2.0)*1.0e-7
    I_stim[u*L_side-1,0] = 0.1
    # Use to stimulate more columns (count from right to left):
#    for u_o in range(1, 5): #1:4
#       I[u*L_side-u_o,0] = 0.1
#       rand_y[u*L_side-1-u_o,1] = (rd.random() - 1.0/2.0)*1.0e-7
#       rand_x[u*L_side-1-u_o,1] = (rd.random() - 1.0/2.0)*1.0e-7

# Iterate time:
for t1 in range(1, t_max+1):
    # Initial conditions:
    if t1 == 1:
        Xo[:,0] = -0.5 + rand_x[:,0]
        Yo[:,0] = -0.5 + rand_y[:,0]
        Zo[:,0] = 0.0
        # Stimulus:
        St = 1 #Count the time steps for stimulation
        I = I_stim
    # Update variables:
    elif t1 > 1:
        Xo = X
        Yo = Y
        Zo = Z
        # Periodic stimulus:
        I = 0.0
        St = St + 1
        if St <= stim_dur:
            I = I_stim
        if St == dt:
            St = 0
    # Diffusive coupling:
    I_snps = diff_coupl(L_side, N_cells, Xo)
    # if t1 == 1:
    #     I_snps = 0.0
    I = I + I_snps
    # Output the membrane potential x(t) and current I(t) from sampling cell:
    time[t1-1,0] = t1
    potential[t1-1,0] = Xo[cell_sample-1,0]
    current[t1-1,0] = I[cell_sample-1,0]
    #Calculate and output the network activity:
    for u in range(0, N_cells): #1:N_cells
        # Activity:
        if Xo[u,0] > 0.0:
            s[t1-1,0] =  s[t1-1,0] + 1
    # Calculate the KTz equations:
    X, Y, Z = ktz_model(N_cells, Xo, Yo, Zo, I)

# Plotting the activity:
plt.plot(time, s, '-o', markersize = 2, linewidth = 1, color = 'blue')
plt.xlabel('t')
plt.ylabel('s(t)')
plt.grid()
plt.show()

# Plot the membrane potential x(t) and current I(t) of one of the central cells:
fig, (ax1, ax2) = plt.subplots(nrows=2, ncols=1, sharex=True,
                               gridspec_kw={'height_ratios': [2.5, 1]})
ax1.plot(time, potential, '-o', markersize = 2, linewidth = 1, color = 'blue')
ax2.plot(time, current, '-o', markersize = 2, linewidth = 1, color = 'red')
ax1.set(ylabel='X(t)')
ax2.set(ylabel='I(t)')
plt.xlabel('t')
plt.xlim(1,t_max)
plt.grid()
plt.show()

# Save network activity to file:
data = np.concatenate((time, s), axis=1)
np.savetxt("activity.dat", data, delimiter=" ", fmt = '%d %.18e')
