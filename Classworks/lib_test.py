import numpy as np
from Library_Classworks import *



A = [[2,1],[1,2]]

B = [[8],[1]]

guess = [[7],[2]]


print("Jacobi Method: ", Gauss_Jacobi_solve(A,B,guess,0.000001))
