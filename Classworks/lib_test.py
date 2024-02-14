from Library_Classworks import *
import numpy as np


def a(x):
    return 4

x,v,t = velocity_verlet_solve(a,0,0,0,1,0.01)

import matplotlib.pyplot as plt
plt.plot(t,x,label='x(t)')
plt.plot(t,v,label='v(t)')
plt.legend()
plt.grid()
plt.show()

