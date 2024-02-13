from Library_Classworks import *
import numpy as np


######################################################################
def f(x):
    return (x**3) - np.cos(x)

# def df(x):
#     return 3*x**2 + np.sin(x)
# print("bisection=",bisection(f,0,1,1e-6,True))
# print("newton_raphson=",newton_raphson_single(f,df,0.5,1e-6))
# print("secant=",secant_method(f,0,1,1e-6))
# print("regula_falsi=",regula_falsi(f,0,1,1e-6,ifbracket=True))
# print("Fixed Point=",fixed_point_single(f,0,1e-6))
######################################################################


######################################################################
# def g1(varlist):
#     return (10 - varlist[0]*varlist[1])**(1/2)
# def g2(varlist):
#     return ((57-varlist[1])/(3*varlist[0]))**(1/2)
# print(fixed_point_multi([g1,g2],[1.5,3.5],1e-6))
######################################################################

# print("Midpoint=",midpoint(f,0,2,100))
# print("Trapezoidal=",trapezoidal(f,0,2,100))
# print("Simpson=",simpsons(f,0,2,100))

######################################################################
# p=Gaussian_Quadrature(f,0,2,6)
# print(p.integrate())
######################################################################


# print("Monte Carlo=",monte_carlo(f,0,2,10000,37))




