"""
Outline:

------One Newton step------
Let n=number of variables.
Suppose x0 has k/2 terms.
Pick random values c1,...,ck.
Solve J(x0(c_i)) * (x1(c_i)-x0(c_i)) = -F(x0(c_i)) for x1(c_i), for all i.
Call these solutions b_i. They are n-tuples.
Do polyfit on {(c_i, b_i[j]) for i=1..k}, for all j, to get the series.
Can do this in one step with polyfit, namely:
pfit(c,transpose([b_1 | b_2 | ... | b_k])).
Output is the nxk matrix of x1, the next result of the Newton step.

Would be nice if it was easy to choose, programmatically, whether we redid the
first k/2 terms, in case this could be a way of picking up on a change in the
initial form as you track the homotopy.

>>> x=np.array([[1],[2],[3]])
>>> x
array([[1],
       [2],
       [3]])
>>> v = np.empty((3,3),int)
>>> v[:,:] = x
>>> v
array([[1, 1, 1],
       [2, 2, 2],
       [3, 3, 3]])
>>> np.multiply.accumulate(v,axis=1)
array([[ 1,  1,  1],
       [ 2,  4,  8],
       [ 3,  9, 27]])

"""
import numpy as np
#from numpy.linalg import lstsq
#from numpy import polyfit as pfit
#------------------------------------------------------------------------------#
def ignoreHeightDot(a,b):
    """ Does 3d matrix mult in horizontal slices """
    return np.dstack([np.matmul(a[i,:,:],b[i,:,:]) for i in xrange(a.shape[0])])
#------------------------------------------------------------------------------#
def newtonStep(x0,F,J,termsTarget,recompute=False):
    """
    Computes next set of terms (coeffs of series).
    x0 is an nxk matrix, n=#variables, k=#terms, a guess for the root.
    F is a list of polynomials.
    J is the Jacobian of F.
    termsTarget is the number of terms to attempt to compute. Becomes k.
        Note: termsTarget shouldn't exceed #terms(x0), probably?
    recompute is whether we re-interpolate for the coeffs of x0 or just reuse
        them.
    """
    k = x0.shape[1] # current number of known terms
    outputSeriesLength = termsTarget + recompute*k
    ci = np.random.rand(termsTarget,1)
    powers = np.empty((outputSeriesLength,outputSeriesLength))
    powers[:,:]=ci
    powers[:,0] = np.array([1 for c in ci])
    np.multiply.accumulate(powers,axis=1,out=powers)
    x0evalled = x0 * powers[:,:k] # componentwise prod, NOT matrix prod
newtonStep(np.array([[1,1,1],[1,1,1]]),0,0,termsTarget=3)

#------------------------------------------------------------------------------#
def alternateNewtonStep(ci,x0evalled,F,J,termsTarget,recompute=False):
    """
    Same as above, except we return x1 evaluated at the ci's plus
    at new points, coupled with the new list of points with our new ones
    appended.
    """
    pass
