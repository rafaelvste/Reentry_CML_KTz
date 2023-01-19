# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt

#KTz model
def ktz_model(X_o, Y_o, Z_o, I_t):
    # Model paramaters
    k = 0.6
    T = 0.154
    Ie = 0.0
    delt = 0.001
    lamb = 0.001
    Xr = -0.48
    # Model equations:
    arg = (X_o - k * Y_o + Z_o + Ie + I_t) / T
    X_f = arg / (1.0 + abs(arg))#Logistic KTz
#    X_f = np.tanh(arg)         #Hyperbolic KTz
    Y_f = X_o
    Z_f = (1.0 - delt) * Z_o - lamb * (X_o - Xr)
    return(X_f, Y_f, Z_f)

# Additional parameters:
t_max = 10000 #Total time steps
dt = 232 #Pacing period
I_stim = 0.1 #Stimulus current value
stim_dur = 10 #Stimulus duration

# Output arrays:
potential = np.zeros((t_max,1))
current = np.zeros((t_max,1))
time = np.zeros((t_max,1), dtype=np.int64)

X = 0.0
Y = 0.0
Z = 0.0

# Iterate time:
for t1 in range(1, t_max+1):
    # Initial conditions:
    if t1 == 1:
        Xo = -0.5
        Yo = -0.5
        Zo = 0.0
        I = I_stim
        St = 1 #Count the timesteps for stimulation
    # Update variables:
    elif t1 > 1:
        Xo = X
        Yo = Y
        Zo = Z
        I = 0.0
        St += 1
        # Periodic stimulus:
        if St <= stim_dur:
            I = I_stim
        if St == dt:
            St = 0
    # Output membrane potential and current:
    time[t1-1,0] = t1
    potential[t1-1,0] = Xo
    current[t1-1,0] = I
    X, Y, Z = ktz_model(Xo, Yo, Zo, I)

#Plotting the membrane potential and current
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

# Save to file
data = np.concatenate((time, potential), axis=1)
np.savetxt("potential.dat", data, delimiter=" ", fmt = '%d %.18e')
