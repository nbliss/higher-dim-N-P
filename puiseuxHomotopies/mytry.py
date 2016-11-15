"""
This is a step of Newton's method on linearized power series.
System: Viviani's curve w/ 2x^2 subbed for x.
Solves (A0 + A1*t)*(B0 + B1*t + ... + B4*t^4) = -F(2,2t)
to give x1-x0.
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
C1 = np.concatenate([A0, Z,Z,Z,Z], axis=1)
C2 = np.concatenate([A1, A0,Z,Z,Z], axis=1)
C3 = np.concatenate([Z, A1,A0,Z,Z], axis=1)
C4 = np.concatenate([Z,Z, A1,A0,Z], axis=1)
C5 = np.concatenate([Z,Z,Z, A1,A0], axis=1)
A = np.concatenate([C1, C2, C3, C4, C5], axis=0)
print A

B = np.array([[a] for a in [0,0,0,0,-4,0,0,0,-4,-4]])
(x10, residuals, rank, singvals) = lstsq(A, B)

reshaped = np.transpose(np.reshape(x10,(5,2)))
x0 = np.array([[0,2,0,0,0],[2,0,0,0,0]])
x1 = x0+reshaped
x1 = x1.tolist()
print "The solution:"
print x1[0]
print x1[1]
print "Residuals:"
print residuals
print "Rank:"
print rank
print "Singular values:"
print singvals
