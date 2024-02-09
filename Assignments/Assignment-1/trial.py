from Library_asgn1 import *

def f(x):
    return np.sqrt(1+x**4)


p=Gaussian_Quadrature(f,0,1,3)
print(p.integrate())