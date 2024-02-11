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
for i in [1]:
    p=Explicit_solve_uxx_ut(g,a,b,2,20,4,5000,i)
    A=p.solve()
    x=np.linspace(0,2,20)
    print(len(A))
    print(A)
    # plt.plot(x,A)
plt.grid()
plt.show()    
