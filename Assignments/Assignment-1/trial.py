from Library_asgn1 import *

######################### Defining the function to be solved #######################
def df(x,y):    # The differential equation is dy/dx = (y+1)^2
    return (y+1)**2
def f(x):       # The analytical solution is y = x/(1-x)
    return x/(1-x)
####################################################################################

#################### Implementing the Fourth Order Runge Kutta #####################
p=ODE_Solve_XY(df,0,0,0.6,1000)
xl,yl=p.RK4_solve()
x_a=np.linspace(0,0.5,1000)
y_a=f(x_a)
import matplotlib.pyplot as plt
plt.plot(xl,yl,label="Fourth Order Runge Kutta",color="red")
plt.plot(x_a,y_a,label="Analytical Solution",color="blue")
plt.legend()
plt.grid()
plt.show()
#####################################################################