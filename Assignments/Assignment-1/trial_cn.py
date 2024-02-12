
import numpy as np
import matplotlib.pyplot as plt
from Library_asgn1 import crank_nicholson

# def crank_nicholson(f, x_0, x_N, N_x, N_t, T, alpha):
#     """"
#     f : function that returns the initial condition
#     x_0 : initial position point
#     x_N : Final position point 
#     N_x : Number of grid points in the x direction
#     N_t : Number of grid points in the t direction
#     T : Final time
#     alpha : Diffusion coefficient
#     """
    
#     # Grid points
#     x = np.linspace(x_0, x_N, N_x+1)
#     t = np.linspace(0, T, N_t+1)
    
#     # Initialize solution matrix
#     u = []
    
#     # Set initial condition (adding to solution matrix)
#     u.append([f(x[i]) for i in range(N_x+1)])
    
#     # set boundary conditions
#     u[0][0] = 0
#     u[0][N_x] = 0
    
#     u = np.array(u)
#     I2paB_inv_cols = [] # list to store the columns of the inverse of (2I+alphaB)
#     for j in range(N_x + 1): # for each col of Identity matrix
#         X = list([0] for i in range(N_x + 1)) # initial guess
#         for step in range(150):
#             flag = 1
#             for i in range(N_x + 1):
#                 sum = 0
#                 if i != 0:
#                     sum += (- alpha * X[i - 1][0])
#                 if i != N_x:
#                     sum += (- alpha * X[i + 1][0])
#                 if i == j:
#                     temp = (1-sum) / (2+(2*alpha))
#                 else:
#                     temp = (-sum) / (2+(2*alpha))
#                 if abs((temp) - (X[i][0])) > tolerance: #checks the tolerance at each new value
#                     flag = 0
#                 X[i][0] = temp
#             if flag == 1:
#                 break
#         if flag == 0:
#             print(f'Eqn not solved after 150 steps for {j}th col')
#             return None
#         I2paB_inv_cols.append(X)
    
#     # building I2paB_inv
#     I2paB_inv = np.array(I2paB_inv_cols[0])
#     for i in range(1,N_x + 1):
#         I2paB_inv = np.append(I2paB_inv,I2paB_inv_cols[i],axis = 1)
    
#     # print(I2paB_inv)
    
#     # multiplication with (2I-alphaB) partially matrix free
#     I2paB_inv_I2naB = []
#     # multiplying with (2I-alphaB)
#     for row in range(N_x + 1):
#         I2paB_inv_I2naB.append([])
#         for col in range(N_x + 1):
#             sum = 0
#             sum += I2paB_inv[row][col] * 2 * (1 - alpha)
#             if col != 0:
#                 sum += I2paB_inv[row][col - 1] * alpha
#             if col != N_x:
#                 sum += I2paB_inv[row][col + 1] * alpha
#             I2paB_inv_I2naB[row].append(sum)
#     I2paB_inv_I2naB = np.array(I2paB_inv_I2naB)
    
#     # print(I2paB_inv_I2naB)    
    
#     del I2paB_inv_cols, I2paB_inv
    
#     # solving the time evolution using Crank-Nicholson method
#     for n in range(N_t):
#         # building the vector (2I-alphaB)u^n
#         u = np.append(u, [I2paB_inv_I2naB @ u[-1,:]], axis = 0) 
#         u[-1,0] = 0
#         u[-1,N_x] = 0

#     ''' 
#     x : Grid points in the x direction
#     t : Grid points in the t direction
#     u : array Solution to the heat equation
#     '''
#     return x,t,u






def u0(x):
    return ((4*x) - ((x**2)/2))

x_0 = 0
x_N = 8
N_x = 80
N_t = 5000
T = 10
alpha1 = (T/N_t)/(((x_N-x_0)/N_x)**2) # alpha = ht/hx^2
alpha2 = 0.5
Xs, Ts, U_xt = crank_nicholson(u0, x_0, x_N, N_x, N_t, T, alpha1)

# plotting the solution
plt.figure(figsize=(10, 6))
plt.plot(Xs, U_xt[0], label=f"t = 0")
for i in range(500, len(Ts), 500):
    plt.plot(Xs, U_xt[i], label=f"t = {Ts[i]:.2f}")
plt.xlabel("x")
plt.ylabel("u(x, t)")
plt.title(f"Solution of the given heat equation using Crank-Nicolson method using alpha = {alpha1:.4f}")
plt.legend()
plt.grid()
plt.show()