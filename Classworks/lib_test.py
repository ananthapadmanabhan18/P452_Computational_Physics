from Library_Classworks import *
import numpy as np


# ilist=[t,y]

def f1(ilist,t):
    return 2*ilist[1]

def f2(ilist,t):
    return ilist[0]

X,Y = shooting_solve([f1,f2],0,1.2,1,0.9,-1.5,0.0001,0.02)

import matplotlib.pyplot as plt
plt.plot(X,Y[1])
plt.show()
