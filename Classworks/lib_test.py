from Library_Classworks import  *

E = 1

def f1(xlist,t):
    return xlist[1]

def f2(ilist,x):
    return ( x**2 - 2*E)*ilist[0]


x,y = shooting_solve([f1,f2],-1,0,1,0,0.0001,0.1,1)



import matplotlib.pyplot as plt
plt.plot(x,y[0])
plt.plot(x,y[1])
plt.show()












