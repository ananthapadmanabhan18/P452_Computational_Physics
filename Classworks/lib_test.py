import numpy as np
from Library_Classworks import *



A = [[2,-1,0],[-1,2,-1],[0,-1,2]]

B = [[-1],[4],[3]]

guess = [[1.2],[0.1],[0.2]]


print(conjugate_gradient(A,B,guess,0.001))

# def f(x):
#     return 4*x**2 + 3*x

# print(simpsons(f,0,10,5000))