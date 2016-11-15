"""
Uses numpy to apply a blocked version of a linear system
solver for a truncated matrix series.
In the example below, we look for a solution of the form
x(t) = x(-1)*t**(-1) + x(0) + x(1)*t.
For the inverse of the matrix, to compute B(-1)*t**(-1) + B(0),
see the script invserse.py, we do not need a square matrix...
"""
import numpy as np
from numpy.linalg import lstsq
A0 = np.array([[0, 4], [0, 0]])
A1 = np.array([[4, 0], [4, 0]])
Z = np.array([[0, 0], [0, 0]])
print('A0 =')
print(A0)
print('A1 =')
print(A1)
C1 = np.concatenate([A0, Z, Z], axis=1)
C2 = np.concatenate([A1, A0, Z], axis=1)
C3 = np.concatenate([Z, A1, A0], axis=1)
A = np.concatenate([C1, C2, C3], axis=0)
print('coefficient matrix A = ')
print(A)
b = np.array([[0], [0], [1], [0], [0], [1]])
print('right hand side b = ')
print(b)
(x, residuals, rank, singvals) = lstsq(A, b)
print('the solution :')
print(x)
print('the singular values :')
print(singvals)
print('the rank : %d' % rank)
print('verification :')
print(np.matrix(A)*np.matrix(x))
