import numpy as np


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


def shooting_solve(fns,x0,y0,x1,y1,tol,h,guess1=0):  
    '''
    # Shooting Method
    ## Parameters
    - fns: List of functions to be solved that is dy_i/dx = f_i(x,t) where x is a list of variables.
    - x0: i nitial value of t or x (acc to qn)
    - y0: v alue of the solution at x0
    - x1: final value of t or x (acc to qn)
    - y1: value of the solution at x1
    - guess1: The initial guess for the second variable
    - tol: The tolerance for the solution
    - h: step size
    ## Returns
    - X: List of x values or t values
    - Y: List of List of y values for each variable y_i
    '''  
    if guess1 == 0:
        guess1 = (y1-y0)/(x1-x0)

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


def finite_element_solve(f: callable,x_i: float, y_i:float,x_f: float,y_f: float,y_f_prime: float, N: int):
    '''
    # Finite Element Method
    for solving the differential equation d2y/dx2 = f(x,y) with boundary conditions y(x_i) = y_i and y'(x_f) = y_f
    ## Parameters
    - f: The function f(x,y) in the differential equation
    - x_i: Initial value of x
    - y_i: Initial value of y
    - x_f: Final value of x
    - y_f: Final value of y
    - y_f_prime: value of y' at x=x_i
    - N: Number of steps to divide the interval [x_i,x_f]
    ## Returns
    - x: List of x values
    - y: List of y values
    '''
    x = np.linspace(x_i, x_f, N+1)
    h = (x_f - x_i) / (N+1)
    A = 2*np.eye(N-1) + np.diagflat([-1 for i in range(N-2)],-1) + np.diagflat([-1 for i in range(N-2)],1)
    A = A.tolist()
    B = []
    for i in range(len(A)):
        if i == 0:
            B.append([(-h**2)*f(x[i],y_i)  + y_i])
        elif i == len(A)-1:
            B.append(f(x[i],y_f) - y_f)
    return A,None
    #not complete








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

def leap_frog_solve(F: callable,pi0: float,x_i: float,x_f: float,t_i: float, t_f: float,N: int):
    '''
    # Leap Frog Method
    Used to solve the hamiltons equation of motion
                    d2x/dt2 = A(x,t)
    '''
    t = np.linspace(t_i, t_f, 2*N)
    dt = (t_f - t_i) / N
    x = np.zeros(2*N)
    pi = np.zeros(2*N)
    x[0] = x_i
    pi[0] = pi0
    pi[1] = pi[0] - dt*F(x[0],t[0])
    x[2] = x[0] + dt*pi[1]

    return pi

# LEAP FROG NOT DONE