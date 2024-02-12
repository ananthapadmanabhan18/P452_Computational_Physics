import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
# Function to solve a tridiagonal system of linear equations using Crank-Nicolson method
def cnicolson(n, alpha):
    # Initialize vector 'v' with values
    v = np.zeros(n)
    for i in range(len(v)):
        v[i] = 4 * i - ((i ** 2) / 2.0)

    # Create identity matrix 'I' and tridiagonal matrix 'B'
    I = np.identity(n)
    B = np.zeros((n, n))
    for i in range(n):
        B[i, i] = 2
    for j in range(n - 1):
        B[j, j + 1] = -1
    for j in range(1, n):
        B[j, j - 1] = -1

    # Construct matrices 'm1' and 'm2' for the Crank-Nicolson method
    m1 = 2 * I - 4 * alpha * B
    m2 = np.linalg.inv(2 * I + 4 * alpha * B)
    print(m2)

    # Initialize vector 'm' with values and an empty list 'l' to store solutions
    m = np.array(v)
    l = []
    n = 0

    # Perform iterations to solve the linear system using Crank-Nicolson method
    while n < 8:
        v = m1 @ m2 @ v
        l.append(v)
        n += 1

    # Return the list of solutions
    return l
# l=cnicolson(20, 0.4)
# print(l)




# Call the cnicolson function with parameters and store the results
stored = cnicolson(8, 0.4)

# Create a Pandas DataFrame from the stored solutions
df = pd.DataFrame(stored)

# Print the Pandas DataFrame to the console
print("Solution as Pandas DataFrame:")
print(df)

# Plot a contour plot using the Pandas DataFrame
plt.figure(figsize=(16, 10))
contour_plot = plt.contourf(df.columns, df.index, df.values, cmap='viridis', levels=20)
plt.colorbar(contour_plot, label='Values')
plt.xlabel('X-axis')
plt.ylabel('Y-axis')
plt.title('Contour Plot from Pandas DataFrame')

# Show the plot
plt.show()