import numpy as np
import scipy as scipy
import math as math
#####################################################################################
#                              Numerical Integration                             
#####################################################################################
################################# Midpoint Method ###################################
def mid_integ(f,a,b,N):                           # f->Function
    h=(b-a)/N                                     # a->Lower Limit
    I=0                                           # b->Upper Limit
    x=a                                           # N->Number of Intervals
    while x<=b:
        I+=h*f(x+h/2)
        x+=h
    return I
################################# Trapezoidal Method ################################
def trap_integ(f,a,b,N):                          # f->Function
    h=(b-a)/N                                     # a->Lower Limit
    I=0                                           # b->Upper Limit
    x=a                                           # N->Number of Intervals
    while x<=b:
        I+=h/2*(f(x)+f(x+h))
        x+=h
    return I
################################# Simpsons Method ###################################
def simp_integ(f,a,b,N):                          # f->Function
    h=(b-a)/N                                     # a->Lower Limit
    I=0                                           # b->Upper Limit
    x=a                                           # N->Number of Intervals
    while x<=b:
        I+=h/3*(f(x)+4*f(x+h/2)+f(x+h))*0.5
        x+=h
    return I
#####################################################################################
#####################################################################################


#####################################################################################
#                            Solution of Linear Equations                             
#####################################################################################
##################################### Bracketing ####################################
def bracket(a0,b0,f):
    n=0
    while f(a0)*f(b0)>=0:
        if abs(f(a0))>abs(f(b0)):
            b0=b0+1.5*(b0-a0)
        else:
            a0=a0-1.5*(b0-a0)       
    return(a0,b0)
################################## Bisection Method #################################
def bijection(a0,b0,f,T):
    a0,b0=bracket(a0,b0,f)
    epsilon=T
    delta=0.001
    count=0
    while (abs(b0-a0))>epsilon:
        c0=(a0+b0)/2
        if f(a0)*f(c0)>0:
            a0=c0
        else:
            b0 = c0 
        count+=1       
    return c0,count 
#################################### Regula Falsi ##################################



########################### Newton-Raphson-Single variable #########################
def fixed_point_single(g,x0,tol):
    x1=g(x0)
    step=1
    while abs(x1-x0)>tol:
        if step>100:
            print("The roots are not converging")
            break
        else:
            x0=x1
            x1=g(x0)
            step+=1
    return x1,step
############################ Newton-Raphson Multi variable ##########################
def fixed_point_multi(glist,x0list,tol):

    if len(glist)!=len(x0list):
        raise IndexError("The number of functions and initial guesses are not equal")
    else:
        for i in range(len(glist)):
            x0list[i] = (glist[i](x0list))
        step=1
        flag=1
        while flag==1:
            if step>100:
                print("The roots are not converging")
                return x0list,step
            else:
                temp = x0list[:]

                for i in range(len(glist)):
                    x0list[i] = (glist[i](x0list))
                step+=1

            for j in range(len(x0list)):
                if abs(temp[j] - x0list[j]) / x0list[j] < tol:
                    flag = 0
                else:
                    flag = 1
                    break
        return x0list,step