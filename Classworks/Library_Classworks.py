import numpy as np
import scipy as scipy
import math as math
from scipy.optimize import root
from Library_old import *
#####################################################################################
#                            Random Number Generation                             
#####################################################################################
class rng():
    def __init__(self,seed, a = 1103515245, c = 12345 ,m = 32768):
        # initiation of data input
        self.term = seed
        self.a = a
        self.c = c
        self.m = m
    def gen(self):
        # generates a random number
        self.term = (((self.a * self.term) + self.c) % self.m)
        return self.term
    def genlist(self,length):
        # returns a list of 'n' random numbers in the range (0,1) where 'n' is 'length'.
        RNs = []
        for i in range(length):
            self.term = (((self.a * self.term) + self.c) % self.m)
            RNs.append(self.term / self.m)
        return RNs  
#####################################################################################
#                        Solution of Non- Linear Equations                             
#####################################################################################
def bracket(a0: float,b0: float,f: float):
    '''
    # Parameters
    - a0: Lower bound of the interval
    - b0: Upper bound of the interval
    # Returns
    - a0: Lower bound of the bracketed interval
    - b0: Upper bound of the bracketed interval
    '''
    n=0
    while f(a0)*f(b0)>=0:
        if abs(f(a0))>abs(f(b0)):
            b0=b0+1.5*(b0-a0)
        else:
            a0=a0-1.5*(b0-a0)       
    return(a0,b0)


def bisection(f: float,a0: float,b0: float,T: float,ifbracket: bool=False):
    '''
    # Bisection Method
    ## Parameters
    - f: Function to find the root     
    - a0: Lower bound of the interval to find the root
    - b0: Upper bound of the interval to find the root
    - T: Tolerance
    - ifbracket: If True, the function will bracket teh interval before using the bisection method
    ## Returns
    - c0: The root of the function
    - count: Number of iterations required tolerance
    '''
    if ifbracket==True:
        a0,b0=bracket(a0,b0,f)
    count=0
    while (abs(b0-a0))>T:
        c0=(a0+b0)/2
        if f(a0)*f(c0)>0:
            a0=c0
        else:
            b0 = c0 
        count+=1      
    return c0,count 



## Newton-Raphson for multivariable is not implemented
def newton_raphson_single(f: float,fd: float,x0: float,T: float):
    '''
    # Newton-Raphson Method for Single Variable
    ## Parameters
    - f: Function to find the root    
    - fd: Derivative of the function f    
    - x0: Initial guess
    - T: Tolerance

    ## Returns
    - xn1: The root of the function
    - count: Number of iterations required to reach the tolerance
    '''
    count=0
    xn=x0
    xn1=xn-(f(xn)/fd(xn))
    while True:
        xn=xn1
        xn1=xn-(f(xn)/fd(xn)) 
        count+=1
        if abs(xn1-xn)<T:
            break          
    return xn1,count+1  


def secant_method(f: float,x0: float,x1: float,tol: float):
    '''
    # Secant Method
    ## Parameters
    - x0: 1st guess
    - x1: 2nd guess
    - f: Function to find the root
    - tol: Tolerance
    ## Returns
    - x2: The root of the function
    - step: Number of iterations required to reach the tolerance
    '''
    x2=x1-((f(x1)*(x1-x0))/(f(x1)-f(x0)))
    step=1
    while abs(x2-x1)>tol:
        if step>100:
            raise ValueError("The roots are not converging")
        else:
            x0=x1
            x1=x2
            x2=x1-f(x1)*(x1-x0)/(f(x1)-f(x0))
            step+=1
    return x2,step

def regula_falsi(f,a0,b0,T,ifbracket=False):
    '''
    # Regula Falsi Method
    ## Parameters
    - f: Function to find the root
    - a0: Lower bound of the interval to find the root
    - b0: Upper bound of the interval to find the root
    - T: Tolerance
    ## Returns
    - c0: The root of the function
    - count: Number of iterations required to reach the tolerance
    '''
    if ifbracket==True:
        a0,b0=bracket(a0,b0,f)
    epsilon=T
    a0,b0=bracket(a0,b0,f)
    for i in range(0,1):
        c0=b0-(((b0-a0)*f(b0))/(f(b0)-f(a0)))
        cold=c0  
        if f(a0)*f(c0)<0:
            b0=c0       
        else:    
            a0=c0 

    count=1

    while True:
        cold=c0
        c0=b0-(((b0-a0)*f(b0))/(f(b0)-f(a0)))
        if f(a0)*f(c0)<0:
            b0=c0       
        else:    
           a0=c0 
        if abs(cold-c0)<epsilon:
            break 
        count+=1
    return c0,count

def fixed_point_single(g: float,x0: float,tol: float):
    '''
    # Fixed Point Method for Single Variable
    ## Parameters
    - g: Function to find the root where the roots to be find is f(x)=0 from where g(x)=x is deducted.
    - x0: Initial guess
    - tol: Tolerance
    ## Returns (x1,step)
    - x1: The root of the function
    - step: Number of iterations required to reach the tolerance
    '''
    x1=g(x0)
    step=1
    while abs(x1-x0)>tol:
        if step>100:
            print("The roots are not converging")
            break
        else:
            x0=x1
            x1=g(x0)
            step+=1
    return x1,step

def fixed_point_multi(glist,x0list,tol):
    '''
    # Fixed Point Method for Multi Variable set of equations
    ## Parameters
    - glist: List of functions to find the root where the roots to be find is f(x)=0 from where g(x)=x is deducted.
    - x0list: List of initial guesses
    - tol: Tolerance
    ## Returns (x0list,step)
    - x0list: The list of roots of the function
    '''
    if len(glist)!=len(x0list):
        raise IndexError("The number of functions and initial guesses are not equal")
    else:
        for i in range(len(glist)):
            x0list[i] = (glist[i](x0list))
        step=1
        flag=1
        while flag==1:
            if step>100:
                print("The roots are not converging")
                return x0list,step
            else:
                temp = x0list[:]

                for i in range(len(glist)):
                    x0list[i] = (glist[i](x0list))
                step+=1

            for j in range(len(x0list)):
                if abs(temp[j] - x0list[j]) / x0list[j] < tol:
                    flag = 0
                else:
                    flag = 1
                    break
        return x0list,step

#####################################################################################
#                              Numerical Integration                             
#####################################################################################
def midpoint(f: float,a: float,b: float,n: int):
    '''
    # Midpoint Method
    ## Parameters
    - f: Function to be integrated
    - a: Lower limit of the integral
    - b: Upper limit of the integral
    - n: Number of intervals
    ## Returns
    - I: The value of the integral
    '''
    h=(b-a)/n
    I=0
    x=a
    while x<=b:
        I+=h*f(x+h/2)
        x+=h
    return I


def trapezoidal(f: float,a: float,b: float,n: int):
    '''
    # Trapezoidal Method
    ## Parameters
    - f: Function to be integrated
    - a: Lower limit of the integral
    - b: Upper limit of the integral
    - n: Number of intervals
    ## Returns (I)
    - I: The value of the integral
    '''
    h=(b-a)/n
    I=0
    x=a
    while x<=b:
        I+=h/2*(f(x)+f(x+h))
        x+=h
    return I

def simpsons(f: float,a: float,b: float,n: int):
    '''
    # Simpson's Method
    ## Parameters
    - f: Function to be integrated
    - a: Lower limit of the integral
    - b: Upper limit of the integral
    - n: Number of intervals
    '''
    h=(b-a)/n
    I=0
    x=a
    while x<=b:
        I+=h/3*(f(x)+4*f(x+h/2)+f(x+h))*0.5
        x+=h
    return I

class Gaussian_Quadrature: 
    '''
    # Gaussian Quadrature (Class)
    ## Parameters
    - f: Function to be integrated
    - a: Lower limit of the integral
    - b: Upper limit of the integral
    - degree: Degree of the LEgendre Polynomial
    ## Returns
    - val: The value of the integral
    '''
    def __init__(self,f,a,b,degree):
        self.a = a
        self.b = b
        self.N = degree
        self.f = f   

    def Pn(self,x,n):
        '''
        # Legendre Polynomial
        ## Parameters
        - x: Variable
        - n: Degree of the polynomial
        ## Returns
        - The value of the Legendre Polynomial
        '''
        if n == 0:
            return 1
        elif n == 1:
            return x
        else:
            return ((2*n-1)*x*self.Pn(x,n-1)-(n-1)*self.Pn(x,n-2))/n

    def Pn_drvtv(self,x,n):
        '''
        # Derivative of Legendre Polynomial
        ## Parameters
        - x: Variable
        - n: Degree of the polynomial
        ## Returns
        - The value of the derivative of the Legendre Polynomial
        '''
        if n == 0:
            return 0
        elif n == 1:
            return 1
        else:
            return (n*(x*self.Pn(x,n)-self.Pn(x,n-1)))/(1-x**2)
    
    def find_legendre_roots(self):
        '''
        # Finding the roots of the Legendre Polynomial
        ## Returns
        - The roots of the Legendre Polynomial of the given degree
        '''
        n=self.N
        num_roots=self.N
        roots = []
        for i in range(1, num_roots + 1):
            guess = np.cos((2*i - 1) * np.pi / (2 * num_roots))
            result = root(self.Pn, guess, args=(n,), jac=self.Pn_drvtv, method='hybr')
            if result.success:
                roots.append(result.x[0])
        return roots
    

    def find_weights(self):
        """
        # Finding the weights for the Gaussian Quadrature   
        ## Returns
        - The weights for the Gaussian Quadrature as a list
        """
        n=self.N
        roots=self.find_legendre_roots()
        weights=[]
        for i in range(n):
            w=2/((1-roots[i]**2)*(self.Pn_drvtv(roots[i],n))**2)
            weights.append(w)
        return weights


    def integrate(self):
        '''
        Algorithm for the Gaussian Quadrature
        '''
        a=self.a
        b=self.b
        n=self.N
        f=self.f
        sum=0
        weights=self.find_weights()
        roots=self.find_legendre_roots()
        for i in range(n):
            y=((b-a)*0.5*roots[i])+((b+a)*.5)
            weightf=weights[i]*f(y)
            sum+=weightf
        val=(b-a)*0.5*sum
        return val     



def monte_carlo(f,a,b,N,seed):
    '''
    # Monte Carlo Integration
    ## Parameters
    - f: Function to be integrated
    - a: Lower limit of the integral
    - b: Upper limit of the integral
    - N: Number of random numbers to be generated
    - seed: Seed for the random number generator
    ## Returns
    - F: The value of the integral
    '''
    p=rng(seed)
    F=0
    for i in range(N):
        k=p.gen()
        k=((b-a)*(k/32768))+a
        F+=((b-a)*f(k))/N   
    return F   

def monte_carlo_error(f,a,b,N,seed):
    '''
    # Monte Carlo Integration
    ## Parameters
    - f: Function to be integrated
    - a: Lower limit of the integral
    - b: Upper limit of the integral
    - N: Number of random numbers to be generated
    - seed: Seed for the random number generator
    '''
    rn=rng(seed)
    F=0
    F1=0
    for i in range(N):
        p=rn.gen()
        p=((b-a)*(p/32768))+a
        F+=f(p)
        F1+=pow(f(p),2)  
    return (F1/N)-pow((F/N),2) 

#####################################################################################
#                              Solving dy/dx = f(x,y)                             
#####################################################################################


