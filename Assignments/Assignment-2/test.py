import numpy as np
from Library_asgn2 import gauss_jordan_solve,LU_Solve_eqn
'''
Defining the Matrix A and B
'''
A=[
    [0,4,2,0,1],
    [4,0,4,10,1],
    [2,5,1,3,13],
    [11,3,0,1,2],
    [3,2,7,1,0]
]
B=[[20],[15],[92],[51],[15]]

'''
Solving using the Gauss-Jordan Method
'''
# Y_GJ_qn2 = gauss_jordan_solve(np.copy(A).tolist(),np.copy(B).tolist())




'''
Solving using the LU Decomposition
'''
Y_LU_qn1 = LU_Solve_eqn(np.copy(A).tolist(),np.copy(B).tolist())

print(np.array(Y_LU_qn1))   

