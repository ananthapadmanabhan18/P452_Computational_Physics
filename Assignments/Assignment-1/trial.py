from Library_asgn1 import *
import numpy as np

N=100



def bracket(a,b,f):
    while f(a)*f(b)>=0:
        if abs(f(a))>abs(f(b)):
            b=b+1.5*(b-a)
        else:
            a=a-1.5*(b-a)       
    return(a,b)

def bisection(f,a,b,tol):
    fn=f
    epsilon=tol
    a,b=Solve_Non_Linear_Equation_1D.bracket(a,b,fn)
    count=0
    while (abs(b-a))>epsilon:
        c=(a+b)/2
        if fn(a)*fn(c)>0:
            a=c
        else:
            b = c 
        count+=1       
    return c,count 






def f(x):
    return (x*np.cos(x)+4*np.sin(x))-((np.cos(n*np.pi/(N+1)))*x)



import matplotlib.pyplot as plt
x= [i for i in range(1,N+1)]
roots=[]
for i in range(1,N+1):
    n=i
    p=Solve_Non_Linear_Equation_1D(f,None,2,2.5,1e-6)
    # print(p.bisection())
    t=p.bisection()
    roots.append(t[0])
    y=[t[0] for i in range(1,N+1)]
    plt.plot(x,y)

plt.show()
print(roots)






# print(fixed_point(f,2.2,1e-6))



