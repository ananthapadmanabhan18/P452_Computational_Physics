from Library_asgn1 import *

def f1(x):
    return x**3-2*x-5



P=Solve_Non_Linear_Equation(f1,2,3,0.000001)
s=P.bisection()
print(s)