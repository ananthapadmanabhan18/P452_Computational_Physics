import numpy as np
import scipy as scipy
import math as math

#####################################################################################
#                            Solving Non-Linear Equations                            
#####################################################################################    
class Solve_Non_Linear_Equation:
    def __init__(self,f,a0,b0,tol):
        self.f = f
        self.a0= a0
        self.b0 = b0
        self.tol = tol

    def bracket(a,b,f):
        while f(a)*f(b)>=0:
            if abs(f(a))>abs(f(b)):
                b=b+1.5*(b-a)
            else:
                a=a-1.5*(b-a)       
        return(a,b)
    
    def bisection(self):
        a=self.a0
        b=self.b0
        fn=self.f
        epsilon=self.tol
        a,b=Solve_Non_Linear_Equation.bracket(a,b,fn)
        count=0
        while (abs(b-a))>epsilon:
            c=(a+b)/2
            if fn(a)*fn(c)>0:
                a=c
            else:
                b = c 
            count+=1       
        return c,count 
    



#####################################################################################
#                              Numerical Integration                             
#####################################################################################
class Numerical_integration:
    def __init__(self, f, a, b,N):
        self.a = a
        self.b = b
        self.N = N
        self.f = f

    def midpoint(self):
        h=(self.b-self.a)/self.N
        I=0
        x=self.a
        while x<=self.b:
            I+=h*self.f(x+h/2)
            x+=h
        return I
    
    def trapezoidal(self):
        h=(self.b-self.a)/self.N
        I=0
        x=self.a
        while x<=self.b:
            I+=h/2*(self.f(x)+self.f(x+h))
            x+=h
        return I


    def simpsons(self):
        h=(self.b-self.a)/self.N
        I=0
        x=self.a
        while x<=self.b:
            I+=h/3*(self.f(x)+4*self.f(x+h/2)+self.f(x+h))*0.5
            x+=h
        return I
    
    def Pn(self,x,n):
        if n == 0:
            return 1
        elif n == 1:
            return x
        else:
            return ((2*n-1)*x*self.Pn(x,n-1)-(n-1)*self.Pn(x,n-2))/n

    def Pn_drvtv(self,x,n):
        if n == 0:
            return 0
        elif n == 1:
            return 1
        else:
            return (n*(x*self.Pn(x,n)-self.Pn(x,n-1)))/(1-x**2)
#####################################################################################
#####################################################################################