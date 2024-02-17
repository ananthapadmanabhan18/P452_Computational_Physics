from Library_Classworks import *
import numpy as np


def g(x):
    return 4*x - (x**2)/2

def a(t):
    return 0

def b(t):
    return 0

import matplotlib.pyplot as plt
ulist=[]
for i in [0,100,200,300,400,500,600,700]:
    x,T,t = crank_nicoson(g,a,b,0,8,0,4,20,5000,i,iflist=True,k=1/4)
    plt.plot(x,T,label='t='+str(i))

plt.legend()
plt.grid()
plt.show()




