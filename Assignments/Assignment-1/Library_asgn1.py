import numpy as np
import scipy as scipy
import math as math
from scipy.optimize import root

#####################################################################################
#                                Matrix Operations                            
#####################################################################################    
class Matrix_Operations:
    '''
    Parameters:
    A: Matrix A
        Matrix is defined as a list of lists
        for example: A = [
                        [1,2,3],
                        [4,5,6],
                        [7,8,9]
                        ]
    '''
    def __init__(self, A):
        self.A = A

    def print_matrix(self,B,flag):
        '''
        This function prints the matrix A defined in the init function
        if flag=1, it prints the matrix B

        '''
        if flag==1:
            for row in B:
                for element in row:
                    print("\t",element, end="\t")
                print()
        elif flag==0:       
            for row in self.A:
                for element in row:
                    print("\t",element, end="\t")
                print()
        else:
            raise ValueError("The flag should be 0 or 1")        


    def matrix_copy(self,A1):
        '''
        This function returns a copy of the matrix A
        '''
        for i in range(len(A1)):
            A1[i]=A1[i][:]
        return A1    



    def add(self,B1):
        '''
        This function adds two matrices A and B
        A is defined in the init function
        B is the parameter of this function
        '''
        A1=self.matrix_copy(self.A)
        if len(A1)!=len(B1) or len(A1[0])!=len(B1[0]):
            raise ValueError("The matrices are not of the same order")
        elif():
            pass
        else:
            for i in range(len(A1)):
                for j in range(len(A1[0])):
                    A1[i][j]=A1[i][j]+B1[i][j]
            return A1
        
    def multiply_by_scalar(self,scalar):
        '''
        The function multiplies a matrix by a scalar
        '''
        A1=self.matrix_copy(self.A)
        for i in range(len(A1)):
            for j in range(len(A1[0])):
                A1[i][j]=A1[i][j]*scalar
        return A1    
        
    def multiply(self,B1):

        '''
        The function multiplies two matrices A and B and the result is BA
        '''
        A1=self.matrix_copy(B1)
        del B1
        B1=self.matrix_copy(self.A)
        if len(A1[0])!=len(B1):
            raise ValueError("Length of the columns of A should be equal to the length of the rows of B")
        else:
            result=[]
            for i in range(len(A1)):
                row=[]
                for j in range(len(B1[0])):
                    sum=0
                    for k in range(len(A1[0])):
                        sum+=A1[i][k]*B1[k][j]
                    row.append(sum)
                result.append(row)
            return result



#####################################################################################
#                            Solving Non-Linear Equations                            
#####################################################################################    
class Solve_Non_Linear_Equation_1D:
    def __init__(self,f,df,a0,b0,tol):
        self.f = f
        self.a0= a0
        self.b0 = b0
        self.tol = tol
        self.df = df
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
    
    def secant(self,guess1,guess2):
        fn=self.f
        t=self.tol
        x0=guess1
        x1=guess2
        x2=x1-((fn(x1)*(x1-x0))/(fn(x1)-fn(x0)))
        step=1
        while abs(x2-x1)>t:
            if step>100:
                raise ValueError("The roots are not converging")
                break
            else:
                x0=x1
                x1=x2
                x2=x1-fn(x1)*(x1-x0)/(fn(x1)-fn(x0))
                step+=1
        return x2,step 
    
    def fixed_point(self,x0):
        g=self.f
        tol=self.tol
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
#####################################################################################
#####################################################################################
                
   
#####################################################################################
#                              Numerical Integration                             
#####################################################################################
class Newton_Cotes:
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

class Gaussian_Quadrature: 
    def __init__(self,f,a,b,degree):
        self.a = a
        self.b = b
        self.N = degree
        self.f = f   

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
    
    def find_legendre_roots(self):
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
        n=self.N
        roots=self.find_legendre_roots()
        weights=[]
        for i in range(n):
            w=2/((1-roots[i]**2)*(self.Pn_drvtv(roots[i],n))**2)
            weights.append(w)
        return weights


    def integrate(self):
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


        
        
        
#####################################################################################
#####################################################################################
    

#####################################################################################
#                                      Solving ODEs                             
#####################################################################################
class ODE_Solve_XY:
    def __init__(self, dy, xi,yi,xf,N):
        self.xi = xi
        self.yi = yi
        self.xf = xf
        self.N = N
        self.dy = dy
        
    def forward_euler(self):
        x_ini=self.yi
        t_ini=self.xi
        t_final=self.xf
        n=self.N
        dx=self.dy
        dt=(t_final-t_ini)/n
        xlist=[]
        ylist=[]
        t=t_ini
        while t<=t_final:
            xlist.append(t)
            ylist.append(x_ini)
            x_ini+=dt*dx(t,x_ini)
            t+=dt
        return xlist,ylist        
    
    def backward_euler(self):
        x0=self.xi
        y0=self.yi
        xf=self.xf
        num_points=self.N
        fn=self.dy

        h = (xf - x0) / num_points
        x_values = np.linspace(x0, xf, num_points + 1)
        y_values = np.zeros(num_points + 1)
        y_values[0] = y0

        for i in range(1, num_points + 1):
            y_values[i] = y_values[i - 1] + h * fn(x_values[i], y_values[i - 1])

        return x_values, y_values

    def predictor_corrector(self):
        dybydx=self.dy
        x0=self.xi
        y0=self.yi
        x_f=self.xf
        n=self.N

        h=(x_f-x0)/n
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
    
    def RK2_solve(self):

        dybydx=self.dy
        x0=self.xi
        y0=self.yi
        xf=self.xf
        N=self.N

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
    
    def RK4_solve(self):

        dybydx=self.dy
        x0=self.xi
        y0=self.yi
        x_f=self.xf
        N=self.N

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

def semi_implicit_euler_solve(f,g,x0,v0,t0,t_max,step_size):
    '''
    Parameters:
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
    Returns:
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

class verlet_algorithm:

    '''
    Parameters:
    - a: a(x)
        The acceleration function
    - x0: initial position: x(t0) = x0
    - v0: initial velocity: v(t0) = v0
    - t0: initial time
    - t_max: final time:
        The to which the Solution is to be calculated
    - step: The size of the interval
    Returns:
    - xlist: List of x values
    - vlist: List of v values
    - tlist: List of t values
    '''

    def __init__(self,a, x0, v0, t0, t_max, step_size):
        self.a = a
        self.x0 = x0
        self.v0 = v0
        self.t0 = t0
        self.t_max = t_max
        self.step_size = step_size

    def verlet_solve(self):
        #Defining the variables
        h=self.step_size
        xlist=[]
        vlist=[]
        tlist=[]
        tm=self.t_max
        x=self.x0
        t=self.t0
        v=self.v0
        acc_fn=self.a

        #The first 
        xlist.append(x)
        vlist.append(v)
        tlist.append(t)
        x1=(x)+(h*v)+(0.5*h*h*self.a(x))
        v1=(x1-x)/h
        t=t+h
        xlist.append(x1)
        vlist.append(v1)
        tlist.append(t)

        #The rest of the steps
        while t<=tm:
            x2=(2*x1)-(x)+(h*h*acc_fn(x1))
            v=(x2-x)/(2*h)
            x=x1
            x1=x2
            t=t+h
            xlist.append(x)
            vlist.append(v)
            tlist.append(t)

        return xlist,vlist,tlist    

    def velocity_verlet(self):
        
        h=self.step_size
        x=self.x0
        v=self.v0
        t=self.t0

        xlist=[x]
        vlist=[v]
        tlist=[t]
        pass

#####################################################################################
#                                      Solving PDEs                             
#####################################################################################
class Explicit_solve_uxx_ut: #not complete
    def __init__(self,g,a,b,L,nx,T,nt,timestep):
        '''
        The equation that we are trying to solve is: u_xx=u_t
        Parameters:
        - g: Function g(x)=u(x,0)
        - a: Function a(x)=u(0,t)
        - b: Function b(x)=u(L,t)
        - L: Length of the rod
        - nx: Number of spatial steps
        - T: Total time
        - nt: Number of time steps
        - timestep: How many times the u(x,t) has to be time-evolved
        '''
        self.g = g
        self.a = a
        self.b = b
        self.L = L
        self.nx = nx
        self.nt = nt
        self.T = T
        self.timestep = timestep
    def check_stability(self):
        '''
        This function checks the stability of the Explicit method
        That is, it checks if the value of alpha is less than 0.5 
        '''
        alpha = (self.T*(self.nx**2))/((self.L**2)*(self.nt))
        if alpha>0.5:
            raise ValueError("The method is not stable")
        else:
            pass

    def create_A_inv(self):
        '''
        This function creates the matrix A
        '''
        self.check_stability()

        alpha = (self.T*(self.nx**2))/((self.L**2)*(self.nt))
        A=[]
        for i in range(self.nx):
            row=[]
            for j in range(self.nx):
                if i==j:
                    row.append(1+2*alpha)
                elif abs(i-j)==1:
                    row.append(-alpha)
                else:
                    row.append(0)
            A.append(row)
        A_inv = np.linalg.inv(A)
        return A_inv 

    def create_V0(self):
        f=self.g
        V0=[]
        for i in range(self.nx):
            if i==0:
                V0.append([self.a(0)])
            elif i==self.nx-1:
                V0.append([self.b(0)])
            else:
                V0.append([f(i*self.L/self.nx)])   
        return V0   

    def solve(self):     
        A=self.create_A()
        V0=self.create_V0()
        A_inv = np.linalg.inv(A)
        for i in range(self.timestep):
            V1=Matrix_Operations(V0).multiply(A_inv)
            # V1=np.matmul(A_inv,V0)
            del V0
            V0=V1
            del V1
            V0[0][0]=0
            V0[self.nx-1][0]=0
        return V0


        



def crank_nicolson(g: float,L: int,n: int, T: float,dt: float):
    """
    - g: g(x) function at t=0 that is g(x) = u(x,0)
    - L: Max value of position
    - n: number of spatial mesh points
    - T: Limit of time
    - dt: time step
    - The initial values a(t) and b(t) are given to be  0 for the case in the question.
    This code is specific for the case in the question.
    """
    dx = L / (n + 1)
    alpha = dt / (dx**2)  
    A = np.zeros((n, n))  
    B = np.zeros((n, n))  
    I = np.zeros((n, n))  

    '''
    The matrix A and B are created as per the Crank Nicolson method,
    the matrix A is a tridiagonal matrix with 2+2*alpha on the diagonal and -alpha on the off-diagonals
    The matrix B is a tridiagonal matrix with 2-2*alpha on the diagonal and alpha on the off-diagonals
    '''
    for i in range(n):
        A[i,i] = 2+2*alpha
        B[i,i] = 2-2*alpha
        for j in range(n):
            if j == i + 1 or j == i - 1:
                A[i,j] = -1*alpha
                B[i,j] = +1*alpha

    x = [0]
    for i in range(n - 1):
        x.append(x[i] + dx)
    '''
    The function g(x) is used to create the 
    initial values of the vector v0'''
    v0 = []
    for i in range(n):
        v0.append(g(x[i]))
    v0[-1] = v0[0]

    v0 = np.array(v0)
    v = v0.copy()
    '''Arrays are initialized to store the solution at each time step'''
    solution_at_each_time = [v0.copy()]
    for _ in range(int(T / dt)):
        C = np.matmul(np.linalg.inv(A), B)
        v = np.matmul(C, v)
        solution_at_each_time.append(v.copy())

    return solution_at_each_time, x



                               

                    



        








        











