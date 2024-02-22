import numpy as np
from Library_Classworks import *

def solve_wave(f0, f_before, xi, xf, ti, tf, nx, nt, alpha = None, zeroboundary = False, times = []):
    # partial differential solver for wave equation.
    # f0 is the intial solution at time = 0
    # xi, xf and ti, tf are the initial x and t respectively.
    # left and right boundary conditions are assumed to be 0
    # nx is the number x divsions and nt is the number of t divisons
    # times is the list of time steps needed on the graphdata
    # returns the graphdata

    hx = (xf-xi)/nx
    ht = (tf-ti)/nt
    if alpha is None:
        alpha = (ht**2)/(hx**2) # alpha = (ht^2)/(hx^2)
        

    # initial V0 data
    V0 = []
    V1 = []
    X = []
    x = xi + hx
    for i in range(1,nx):
        X.append(x)
        V0.append(f0(x))
        V1.append(f0(x))
        x += hx
    Vlist = [V0[:]]
    
    if f_before is None or zeroboundary is True: # assume zero boundary
        V_before = [V0[0] + ((alpha / 2) * (V0[1] - (2 * V0[0])))]
        for i in range(2, nx-1):
            V_before.append(V0[i] + ((alpha / 2) * (V0[i+1] - (2 * V0[i]) + V0[i-1])))
        V_before.append(V0[nx-1] + ((alpha / 2) * (-(2 * V0[nx-1]) + V0[nx-2])))
    else:
        V_before = [f_before(x) for x in X]
    
    for i in range(1,nt+1):
        V1 = [2*V0[0] - V_before[0] + (alpha * (V0[1] - (2 * V0[0])))]
        for j in range(1,nx-2):
            V1.append(2*V0[j] - V_before[j] + (alpha * (V0[j + 1] - (2 * V0[0]) + V0[j - 1])))
        V1.append(2*V0[nx-2] - V_before[nx-2] + (alpha * (-(2 * V0[nx-2]) + V0[nx-3])))
        
        V_before = V0[:]
        V0 = V1[:]

        if i in times:
            Vlist.append(V1[:])

    V0.insert(0,0)
    V0.append(0)
    X.insert(0,xi)
    X.append(x) #currently x is the last value, xf
    
    if bool(times) is False:
        return X,V0
    
    for i in range(len(Vlist)):
        Vlist[i].insert(0,0)
        Vlist[i].append(0)
    
    return X,Vlist



def f0(x):
    return 2*np.sin(np.pi*x)

def f_before(x):
    return 2*np.sin(np.pi*x)

xi = 0
xf = 4
ti = 0
tf = 2
nx = 100
nt = 10
alpha = None
zeroboundary = False

X,Vlist = solve_wave(f0, f_before, xi, xf, ti, tf, nx, nt, alpha, zeroboundary,[0,0.1,0.2])

# print(Vlist)

import matplotlib.pyplot as plt

# plt.plot(X,Vlist[0])

for i in [0.1,0.2,0.3]:
    X,Vlist = solve_wave(f0, f_before, xi, xf, ti, i, nx, nt, alpha, zeroboundary,[0])
    plt.plot(X,Vlist[0])

plt.show()    




    
