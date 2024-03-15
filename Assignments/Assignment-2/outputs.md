<h1> Question - 1 Outputs</h1>


A is a symmetric Matrix,So we can proceed with Cholesky Decomposition and Gauss Seidel Method
The solution of the syatem of the linear equations are tabluated below:

+----------+----------+--------------+
| Variable | Cholesky | Gauss Seidel |
+----------+----------+--------------+
|   x_1    | 1.000000 |  0.9999998   |
|   x_2    | 1.000000 |  0.9999998   |
|   x_3    | 1.000000 |  0.9999999   |
|   x_4    | 1.000000 |  0.9999999   |
|   x_5    | 1.000000 |  0.9999999   |
|   x_6    | 1.000000 |  0.9999999   |
+----------+----------+--------------+
The Number of Iterations for Gauss Seidel Method is:  16


<h1> Question - 2 Outputs</h1>

The solution of the syatem of the linear equations are tabluated below:

+----------+--------------+------------------+
| Variable | Gauss Jordan | LU Decomposition |
+----------+--------------+------------------+
|   x_1    |    2.6746    |      2.6746      |
|   x_2    |    3.7119    |      3.7119      |
|   x_3    |   -0.0533    |     -0.0533      |
|   x_4    |   -0.0744    |     -0.0744      |
|   x_5    |    5.2591    |      5.2591      |
+----------+--------------+------------------+



<h1> Question - 3 Outputs</h1>

The Number of Iterations for Conjugate Gradient Method is 6 which is same as the length of the A
The solution of the syatem of the linear equations are tabluated below:

+----------+----------+
| Variable | Solution |
+----------+----------+
|   x_1    | -0.7245  |
|   x_2    |  0.2177  |
|   x_3    |  0.9666  |
|   x_4    | -0.6857  |
|   x_5    | -0.0382  |
|   x_6    |  0.6488  |
+----------+----------+




The matrix A is shown below:

[[ 2. -1.  0.  0.  0.  0.]
 [-1.  4. -1.  0. -1.  0.]
 [ 0. -1.  4.  0.  0. -1.]
 [ 0.  0.  0.  2. -1.  0.]
 [ 0. -1.  0. -1.  4. -1.]
 [ 0.  0. -1.  0. -1.  4.]]


The Number of Iterations for Conjugate Gradient Method is 6 which is same as the length of the A
The inverse of the Matrix A is shown below:

[[0.5868 0.1735 0.0501 0.0286 0.0572 0.0268]
 [0.1735 0.347  0.1002 0.0572 0.1145 0.0537]
 [0.0501 0.1002 0.297  0.0268 0.0537 0.0877]
 [0.0286 0.0572 0.0268 0.5868 0.1735 0.0501]
 [0.0572 0.1145 0.0537 0.1735 0.347  0.1002]
 [0.0268 0.0537 0.0877 0.0501 0.1002 0.297 ]]


For the verification of the result the product of A and A_inv is shown below:

[[ 1.  0. -0. -0.  0. -0.]
 [-0.  1.  0. -0. -0. -0.]
 [-0. -0.  1.  0.  0.  0.]
 [-0.  0. -0.  1.  0. -0.]
 [-0. -0. -0. -0.  1.  0.]
 [-0.  0. -0.  0.  0.  1.]]