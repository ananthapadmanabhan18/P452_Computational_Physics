import numpy as np
from Library_Classworks import *



def F(x):
    return -x

x_i = 0
x_f = 10
t_i = 0
t_f = 10
N = 100


t,x,pi = leap_frog_solve(F, 4,x_i, x_f, t_i, t_f, N)


import matplotlib.pyplot as plt
plt.plot(t,x)
plt.show()