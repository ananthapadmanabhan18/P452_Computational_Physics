from Library_Classworks import *
import numpy as np


def g(x):
    return 20*abs(np.sin(np.pi*x))

def a(t):
    return 0

def b(t):
    return 0

x,T,t = pde_explicit_solve(g,a,b,0,2,0,4,20,5000,100)
# import matplotlib.pyplot as plt

# plt.contourf(x, t, T)
# plt.colorbar()
# plt.xlabel('x')
# plt.ylabel('t')
# plt.title('Contour Plot of T')
# plt.show()
