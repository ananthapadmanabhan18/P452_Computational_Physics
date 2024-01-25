import numpy as np
import scipy as sp


def pde_explicit(f,nx,nt,lx,lt,N):
    hx=lx/nx
    ht=lt/nt
    a = ht/(pow(hx,2))
    V0,V1 = [0],[0]
    for i in range(nx+1):
        V1.append(0)
    for i in range(nx+1):
        V0.append(f(i*hx))
    for j in range(N):
        for i in range(nx+1):
            if i==0:
                V1[i]=(1-2*a)*V0[i] + (a*V0[i+1])
            elif i==nx:
                V1[i]=(a*V0[i-1]) + (1-2*a)*V0[i]
            else:
                V1[i]=(a*V0[i-1]) + (1-2*a)*V0[i] + (a*V0[i+1])
        for i in range(nx+1):
            V0[i] = V1[i]        
    if N==0:
        return V0
    else:
        return V1




def gen_RK4_coupled(fnlist,x0,y0s,limit,h):
    limit -= h/2 
    n = len(y0s) 
    k1 = [0 for i in range(n)]
    k2 = [0 for i in range(n)]
    k3 = [0 for i in range(n)]
    k4 = [0 for i in range(n)]
    tys= [0 for i in range(n)] 
    datY = []
    for i in range(n):
        datY.append([y0s[i]])
    datT = [x0]
    while x0 < limit:
        for i in range(n):
            k1[i] = h*fnlist[i](y0s,x0)    
        for i in range(n):
            tys[i] = y0s[i] + (k1[i] / 2)
        for i in range(n):
            k2[i] = h*fnlist[i](tys, (x0 + (h/2)))
        for i in range(n):
            tys[i] = y0s[i] + (k2[i] / 2)
        for i in range(n):
            k3[i] = h*fnlist[i](tys, (x0 + (h/2)))   
        for i in range(n):
            tys[i] = y0s[i] + k3[i]
        for i in range(n):
            k4[i] = h*fnlist[i](tys, (x0 + h))
        for i in range(n):
            y0s[i] += ((k1[i] + (2 * k2[i]) + (2 * k3[i]) + (k4[i])) / 6)
        x0 += h
        for i in range(n):
            datY[i].append(y0s[i])
        datT.append(x0)
    return datT, datY