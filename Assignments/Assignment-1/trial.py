from Library_asgn1 import *

def a(x):
    return -x

p=verlet_algorithm(a,0,3,0,5,0.01)
x,v,t = p.verlet_solve()

import matplotlib.pyplot as plt
plt.plot(t,x,label="x")
plt.plot(t,v,label= "v")
plt.show()