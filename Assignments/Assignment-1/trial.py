from Library_asgn1 import *
import numpy as np



def f(x):
    return 20*abs(np.sin(np.pi*x))

w=Explicit_solve_uxx_ut(f,0,0,2,100,20,int(50000))
bleh=[]
for i in range(len(w)):
    bleh.append(w[i][0])


import matplotlib.pyplot as plt
x=np.linspace(0,2,20)
plt.plot(x,bleh)
plt.show()
