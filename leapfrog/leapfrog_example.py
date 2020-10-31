# !/usr/bin/python3
#  -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
import numpy as np


# Constants and canonical units
Pi = np.pi
Mcan, Rcan, G = 1.989E33, 6.957E10, 1
Tcan = np.sqrt(Rcan**3/(6.674E-8*Mcan))                   
Vcan = Rcan/Tcan


# Calculates magnitude of vector=[v1,v2,...]
def mag(vector):
    return np.sqrt(vector.dot(vector))


# Calculates direction of forces
def direc(X_e,X_s):
    return np.array([(X_s[0]-X_e[0])/mag(X_e-X_s), (X_s[1]-X_e[1])/mag(X_e-X_s)])


# Aceleration function (for the main proyect this will be mooooore complicated)
# Acc = [A_e, A_s] = [[A1x, A1y, A1z], [A2x, A2y, A2z]]
def grav(r,M,H):
    if r < np.sum(H):
        Acc = np.ones(2)*M[::-1]/np.sum(H)**2
    else:
        Acc = np.ones(2)*M[::-1]/r**2
    return Acc


# Timestep calculates 0.25*dt*sqrt(size of particle/mag(Aceleration)) for each particle
# and spits out the minimum of these and dt.
def timestep(dt,A,H):
    dta = 0.25*dt*min([np.sqrt(H[i[0]]/mag(i[1])) for i in enumerate(A)])
    return min(dta,dt)

def leapfrog(dt,X,V,A):
# Input form:
# X = [X_e, X_s] ===> shape: (2,2)  (Number of particles, number of components)
# V = [V_e, V_s] ===> shape: (2,2)  (Number of particles, number of components)
# A = [A_e, A_s] ===> shape: (2,2)  (Number of particles, number of components)


# Position integration 
    X = X + V*dt + 0.5*A*dt**2

# Aceleration update
    Aupdate = grav(mag(X[0]-X[1]), [M_e,M_s], [h_e,h_s])
    Aupt = np.ones((2,2))
    Aupt[0] = Aupdate[0]*direc(X[0], X[1])                                  # Components of Earth aceleration
    Aupt[1] = -Aupdate[1]*direc(X[0], X[1])                                 # Components of Sun aceleration
    
# Velocity integration    
    V = V + 0.5*(A + Aupt)*dt
    return [X, V, Aupt]

X, V, T = [[],[],[],[]], [[],[],[],[]], []

# Masses and sizes of Sun and Earth
M_s, M_e = 1, 3E-6
h_s, h_e = 1, 9.168E-3


# Initial conditions for Sun-Earth orbit
X0_e, V0_e = np.array([1.521E13/Rcan, 0]), np.array([0, 2.93E6/Vcan])
X0_s, V0_s = np.array([0, 0]), np.array([0, 0])

X0_e, V0_e = np.array([1.521E13/Rcan, 0]), np.array([-2e5/Vcan, 9.3E5/Vcan])

dt0 = 0.2                                                                   # This will be inside the loop for dt variable
t = 0
[A0_e, A0_s] = grav(mag(X0_e-X0_s), [M_e,M_s], [h_e,h_s])
X_e, X_s, V_e, V_s = X0_e, X0_s, V0_e, V0_s                                 # Setting initial conditions for integration loop
A_e = A0_e*direc(X_e, X_s)                                                  # Components of Earth aceleration        
A_s = -A0_s*direc(X_e, X_s)                                                 # Components of Sun aceleration

ndump = 0

Xt, Vt = np.array([X_e, X_s]), np.array([V_e, V_s])                         # Order the inital conditions in the used format 
At = np.array([A_e, A_s])                                                   #

# Integration loop
while t < 1.98E4:
    if t >= ndump:
        ndump+=10                                                           # Stores T,X,V,A every ndumps (time units)
        for i in enumerate(np.hstack(Xt)):                                  #
           X[i[0]] = np.append(X[i[0]],i[1])                                #
        for i in enumerate(np.hstack(Vt)):                                  #
           V[i[0]] = np.append(V[i[0]],i[1])                                #
        T.append(t)                                                         #
        print(t,dt0,timestep(dt0,At,[h_e, h_s]))
    dt = timestep(dt0,At,[h_e, h_s])
    # dt = dt0
    [Xt, Vt, At] = leapfrog(dt, Xt, Vt, At) 
    t+=dt

K = 0.5*(M_e*(V[0]**2 + V[1]**2) + M_s*(V[2]**2 + V[3]**2))                 # Kinetic energy from stored data
U = -M_e*M_s/np.sqrt((X[3]-X[1])**2 + (X[2]-X[0])**2)                       # Potential energy from stored data

Err = np.abs((K+U-K[0]-U[0])/(U[0]+K[0]))                                   # Relative error for total energy 


plt.plot(X[0],X[1],'g',ms=3)
plt.plot(X[2],X[3],'orange',ms=3)
plt.figure(2)
plt.semilogy(T,Err)
plt.show()


