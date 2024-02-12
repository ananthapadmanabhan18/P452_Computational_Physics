from Library_asgn1 import *
import numpy as np



def g(x):
    return 20*abs(np.sin(np.pi*x))
def a(t):
    return 0
def b(t):
    return 0

B=[]
import matplotlib.pyplot as plt
for i in [0,10,20,40]:
    p=crank_nicolson(g,a,b,2,20,4,5000,i)
    hx=2/21
    x,A=p.solve()
    plt.plot(x,A)
print(x)    
plt.grid()
plt.show()    
