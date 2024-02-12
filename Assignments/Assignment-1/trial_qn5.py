from Library_asgn1 import *
import numpy as np
import pandas as pd


def f(x):
    return (4*x) - (x**2)/2


ans = crank_nicolson(f,8,50,4,0.1)

ulist,x=ans


import matplotlib.pyplot as plt
num_times = len(ulist)
num_x_values = len(x)

X, Y = np.meshgrid(x, range(num_times))
u_values = np.array(ulist)
plt.figure(figsize=(8, 6))
contour = plt.contourf(X, Y, u_values, cmap='viridis')
plt.colorbar(contour, label='u')
plt.xlabel('x')
plt.ylabel('Time')
plt.title('Contour Plot of u Values')
plt.show()

