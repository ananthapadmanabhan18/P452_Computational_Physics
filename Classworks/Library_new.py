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



