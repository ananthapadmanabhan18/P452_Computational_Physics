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



def monte_carlo(f: float,a: float,b: float,N: int,seed: int):
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

def monte_carlo_error(f: float,a: float,b: float,N: int,seed: int):
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
def forward_euler(dx: float,x_ini: float,t_ini: float,t_final: float,N: int):
    '''
    # Forward Euler Method
    for solving the differential equation dy/dx = f(x,y)
    also called explicit Euler method
    ## Parameters
    - dx: The function f(x,y): dy/dx = f(x,y)
    - x_ini: Initial value of y such that y(t_ini) = x_ini
    - t_ini: Initial value of x
    - t_final: Final value of x
    - N: Number of steps to divide the interval [t_ini,t_final]
    ## Returns
    - xlist: List of x values
    - ylist: List of y values satisfying the function dy/dx = f(x,y)
    '''
    dt=(t_final-t_ini)/N
    xlist=[]
    ylist=[]
    t=t_ini
    while t<=t_final:
        xlist.append(t)
        ylist.append(x_ini)
        x_ini+=dt*dx(t,x_ini)
        t+=dt
    return xlist,ylist

def backward_euler(f: float,y0: float,x0: float,xf: float,num_points: int):
    '''
    # Backward Euler Method
    for solving the differential equation dy/dx = f(x,y)
    ## Parameters
    - f: The function f(x,y): dy/dx = f(x,y)
    - y0: Initial value of y such that y(x0) = y0
    - x0: Initial value of x
    - xf: Final value of x
    - num_points: Number of steps to divide the interval [x0,xf]
    ## Returns
    - x_values: List of x values
    - y_values: List of y values satisfying the function dy/dx = f(x,y)
    '''
    h = (xf - x0) / num_points
    x_values = np.linspace(x0, xf, num_points + 1)
    y_values = np.zeros(num_points + 1)
    y_values[0] = y0

    for i in range(1, num_points + 1):
        # Use backward Euler method formula: y[i] = y[i-1] + h * f(x[i], y[i])
        y_values[i] = y_values[i - 1] + h * f(x_values[i], y_values[i - 1])

    return x_values, y_values



def predictor_corrector(dybydx: float,y0: float,x0: float,x_f: float,N: int):
    '''
    # Predictor Corrector Method
    for solving the differential equation dy/dx = f(x,y)
    ## Parameters
    - dybydx: The function f(x,y): dy/dx = f(x,y)   
    - y0: Initial value of y such that y(x0) = y0
    - x0: Initial value of x
    - x_f: Final value of x
    - N: Number of steps to divide the interval [x0,xf]
    ## Returns
    - xlist: List of x values
    - ylist: List of y values satisfying the function dy/dx = f(x,y)
    '''
    h=(x_f-x0)/N
    xlist=[]
    ylist=[]
    x=x0
    y=y0
    xlist.append(x)
    ylist.append(y)
    while x<x_f:
        k1=dybydx(x,y)*h
        k2=dybydx(x+h,y+k1)*h
        y=y+0.5*(k1+k2)
        x=x+h
        xlist.append(x)
        ylist.append(y)
    return xlist,ylist

def RK2_solve(dybydx: float,y0: float,x0: float,xf: float,N: int):
    '''
    # Runge-Kutta 2nd Order Method
    for solving the differential equation dy/dx = f(x,y)
    ## Parameters
    - dybydx: The function f(x,y): dy/dx = f(x,y)
    - y0: Initial value of y such that y(x0) = y0
    - x0: Initial value of x
    - xf: Final value of x
    - N: Number of steps to divide the interval [x0,xf]
    ## Returns
    - xlist: List of x values
    - ylist: List of y values satisfying the function dy/dx = f(x,y)
    '''
    h=(xf-x0)/N
    xlist=[]
    ylist=[]
    x=x0
    y=y0
    xlist.append(x)
    ylist.append(y)
    while x<xf:
        k1=h*dybydx(x,y)
        k2=dybydx(x+(h/2),y+(k1/2))*h
        y=y+k2
        x=x+h
        xlist.append(x)
        ylist.append(y)
    return xlist,ylist

def RK4_solve(dybydx: float,y0: float,x0: float,x_f: float,N: int):
    '''
    # Runge-Kutta 4th Order Method
    for solving the differential equation dy/dx = f(x,y)
    ## Parameters
    - dybydx: The function f(x,y): dy/dx = f(x,y)
    - y0: Initial value of y such that y(x0) = y0
    - x0: Initial value of x
    - x_f: Final value of x
    - N: Number of steps to divide the interval [x0,xf]
    ## Returns
    - xlist: List of x values
    - ylist: List of y values satisfying the function dy/dx = f(x,y)
    '''
    h=(x_f-x0)/N
    xlist=[]
    ylist=[]
    x=x0
    y=y0
    xlist.append(x)
    ylist.append(y)
    while x<x_f:
        k1=h*dybydx(x,y)
        k2=h*dybydx(x+(h/2),y+(k1/2))
        k3=h*dybydx(x+(h/2),y+(k2/2))
        k4=h*dybydx(x+h,y+k3)
        y=y+(k1+2*k2+2*k3+k4)/6
        x=x+h
        xlist.append(x)
        ylist.append(y)
    return xlist,ylist 

def RK4_solve_coupled(fnlist,x0,y0s,limit,h):
    '''
    # Runge-Kutta 4th Order Method for Coupled Equations
    ## Parameters
    - fnlist: List of functions to be solved that is dy_i/dx = f_i(x,t) where x is a list of variables.
    - x0: the initial value of t or x (acc to qn)
    - y0s: the value of each of the solution at x0
    - limit: The limit to which the plot is to be made
    - h: step size
    ## Returns
    - datT: List of x values or t values
    - datY: List of List of y values for each variable y_i
    -
    '''
    limit -= h/2 
    n = len(y0s) 
    k1 = [0 for i in range(n)]
    k2 = [0 for i in range(n)]
    k3 = [0 for i in range(n)]
    k4 = [0 for i in range(n)]
    tys= [0 for i in range(n)] 
    datY = []
    for i in range(n):
        datY.append([y0s[i]])
    datT = [x0]
    while x0 < limit:
        for i in range(n):
            k1[i] = h*fnlist[i](y0s,x0)    
        for i in range(n):
            tys[i] = y0s[i] + (k1[i] / 2)
        for i in range(n):
            k2[i] = h*fnlist[i](tys, (x0 + (h/2)))
        for i in range(n):
            tys[i] = y0s[i] + (k2[i] / 2)
        for i in range(n):
            k3[i] = h*fnlist[i](tys, (x0 + (h/2)))   
        for i in range(n):
            tys[i] = y0s[i] + k3[i]
        for i in range(n):
            k4[i] = h*fnlist[i](tys, (x0 + h))
        for i in range(n):
            y0s[i] += ((k1[i] + (2 * k2[i]) + (2 * k3[i]) + (k4[i])) / 6)
        x0 += h
        for i in range(n):
            datY[i].append(y0s[i])
        datT.append(x0)
    return datT, datY


def shooting_solve(fns,x0,y0,x1,y1,guess1,tol,h):  
    '''
    # Shooting Method
    ## Parameters
    - fns: List of functions to be solved that is dy_i/dx = f_i(x,t) where x is a list of variables.
    - x0: first initial value of t or x (acc to qn)
    - y0: first value of the solution at x0
    - x1: second final value of t or x (acc to qn)
    - y1: second value of the solution at x1
    - guess1: The initial guess for the second variable
    - tol: The tolerance for the solution
    - h: step size
    ## Returns
    - X: List of x values or t values
    - Y: List of List of y values for each variable y_i
    '''  
    X,Y = RK4_solve_coupled(fns,x0,[y0,guess1],x1,h)
    ye1 = Y[0][-1]
    if abs(ye1 - y1) < tol:
        return X,Y
    if ye1 < y1:
        guess1side = -1
    else :
        guess1side = 1
    guess2 = guess1 + 2   
    X,Y =RK4_solve_coupled(fns,x0,[y0,guess2],x1,h)
    ye2 = Y[0][-1]
    if ye2 < y1:
        guess2side = -1
    else :
        guess2side = 1
    while guess1side * guess2side != -1:

        if abs(y1-ye2) > abs(y1-ye1):
            guess2 = guess1 - abs(guess2-guess1)
        else:
            guess2 += abs(guess2-guess1)
        X,Y = RK4_solve_coupled(fns,x0,[y0,guess2],x1,h)
        ye2 = Y[0][-1]
        if ye2 < y1:
            guess2side = -1
        else :
            guess2side = 1
    i = 0
    while True:
        newguess = guess1 + (((guess2 - guess1)/(ye1 - ye2))*(y1 - ye2))
        i += 1
        X,Y = RK4_solve_coupled(fns,x0,[y0,newguess],x1,h)
        yvalnew = Y[0][-1]
        if abs(yvalnew - y1) <tol:
            break
        if yvalnew < y1:
            guess1 = newguess
            ye1 = yvalnew
        else:
            guess2 = newguess
            ye2 = yvalnew
    return X,Y



def semi_implicit_euler_solve(f,g,x0,v0,t0,t_max,step_size):
    '''
    # Semi-Implicit Euler Method
    ## Parameters:
    - f: f(v,t):
        dx/dt = f(v,t)
    - g: g(x,t)
        dv/dt = g(x,t)
    - x0: initial position: x(t0) = x0
    - v0: initial velocity: v(t0) = v0
    - t0: initial time
    - t_max: final time:
        The to which the Solution is to be calculated
    - step: The size of the interval
    ## Returns:
    - xlist: List of x values
    - vlist: List of v values
    - tlist: List of t values
    The function returns a 3-tuple of lists containing the x, v and t values respectively
    '''
    h=step_size
    vlist=[]
    xlist=[]
    tlist=[]
    x=x0
    v=v0
    t=t0
    while t<=t_max:
        xlist.append(x)
        vlist.append(v)
        tlist.append(t)
        v=v + (h*g(x,t))
        x=x + (h*f(v,t))
        t=t+h
    return xlist,vlist,tlist    


def verlet_solve(a: float,x0: float,v0: float,t0: float,t_max: float,h: float):
    '''
    # Verlet Method
    ## Parameters
    - a: The function a(x) which gives the acceleration of a(x(t)) = F(x(t))/m
    - x0: initial position: x(t0) = x0
    - v0: initial velocity: v(t0) = v0
    - t0: initial time
    - t_max: final time:
        The to which the Solution is to be calculated
    - h: The size of the interval that (t0,t_max) is divided into
    ## Returns   
    - xlist: List of x values
    - vlist: List of v values
    - tlist: List of t values 
    '''
    xlist=[]
    vlist=[]
    tlist=[]
    x=x0
    t=t0
    v=v0
    xlist.append(x)
    vlist.append(v)
    tlist.append(t)
    x1=(x)+(h*v)+(0.5*h*h*a(x))
    v1=(x1-x)/h
    t=t+h
    xlist.append(x1)
    vlist.append(v1)
    tlist.append(t)
    while t<=t_max:
        x2=(2*x1)-(x)+(h*h*a(x1))
        v=(x2-x)/(2*h)
        x=x1
        x1=x2
        t=t+h
        xlist.append(x)
        vlist.append(v)
        tlist.append(t)
    return xlist,vlist,tlist    

def velocity_verlet_solve(a: callable,x0: float,v0: float,t0: float,t_f: float,h: float):
    '''
    # Velocity Verlet Method
    ## Parameters
    - a: The function a(x) which gives the acceleration of a(x(t)) = F(x(t))/m
    - x0: initial position: x(t0) = x0
    - v0: initial velocity: v(t0) = v0
    - t0: initial time
    - t_f: final time:
        The to which the Solution is to be calculated
    - h: The size of the interval that (t0,t_max) is divided into
    ## Returns
    - xlist: List of x values
    - vlist: List of v values
    - tlist: List of t values    
    '''
    xlist=[x0]
    vlist=[v0]
    tlist=[t0]
    i=1
    while t0<t_f:
        v_half=v0+(0.5*h*a(x0))
        x1= x0 +v_half*h
        v1= v_half+(0.5*h*a(x1))
        xlist.append(x1)
        vlist.append(v1)
        t0+=h
        tlist.append(t0)
        x0=x1
        v0=v1
        i+=1
    return xlist,vlist,tlist 

def leap_frog_solve(pi: callable,F: callable,h):
    '''
    # Leap Frog Method
    Used to solve the hamiltons equation of motion
    '''
    pass   
# LEAP FROG NOT DONE











