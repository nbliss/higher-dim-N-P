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
↑↑↑↑ Not exactly - we want to do a 3d version of J*C*x = F


Would be nice if it was easy to choose, programmatically, whether we redid the
first k/2 terms, in case this could be a way of picking up on a change in the
initial form as you track the homotopy.
"""
import numpy as np
#from numpy.linalg import lstsq
#from numpy import polyfit as pfit
#------------------------------------------------------------------------------#
def ignoreHeightDot(a,b):
    """ Does 3d matrix mult in horizontal slices """
    return np.dstack([np.matmul(a[i,:,:],b[i,:,:]) for i in xrange(a.shape[0])])
#------------------------------------------------------------------------------#
def newtonStep(x0,F,J,termsTarget):
    """
    Computes next set of terms (coeffs of series).
    x0 is an nxk matrix, n=#variables, k=#terms, a guess for the root.
    F is a list of polynomials.
    J is the Jacobian of F.
    termsTarget is the number of terms to attempt to compute. Becomes k.
        Note: termsTarget shouldn't exceed 2 * #terms(x0), probably?
    For this one, we're just going to recompute the first k terms.
    """
    R = F[0].parent()
    print "x0:\n",x0
    k = x0.shape[1] # current number of known terms
    #ci = np.random.rand(termsTarget,1)
    # fixed seed for predictability:
    ci = np.array([[i*0.1] for i in xrange(termsTarget)])
    ci = np.array([[i] for i in xrange(2,termsTarget+2)])
    print "c_i's:\n",ci
    # powers will be rows of ci^column#, i.e. rows are ci, columns are term#
    powers = np.empty((termsTarget,termsTarget))
    powers[:,:]=ci
    powers[:,0] = np.array([1 for c in ci])
    np.multiply.accumulate(powers,axis=1,out=powers)
    print "Powers:\n",powers
    # for each row in x0, take dot products with rows of powers:
    # rows are ci's, columns are different series components (variables)
    x0evalled = np.dot(powers[:,:k],np.transpose(x0))
    Jevalled = []
    Fevalled = []
    for i in xrange(termsTarget):
        #subDict = {R.gens()[j+1]: x0evalled[i][j] for j in xrange(R.ngens()-1)}
        subDict = {xi: x0ci for xi,x0ci in zip(R.gens()[1:],x0evalled[i])}
        subDict[R.0] = ci[i][0] # for first variable, sub just the value
        Ji = np.array(map(lambda a:a.subs(subDict),J.list()))
        Ji.shape = (len(subDict)-1,len(subDict)-1)
        Fi = np.array(map(lambda f:[f.subs(subDict)],F)) #col vector of F@x0@ci
        Jevalled.append(Ji)
        Fevalled.append(Fi)
    J = np.dstack(Jevalled) #eqs x #vars x #ci's
    F = np.dstack(Fevalled).transpose() #cis x 1 x #vars
    C = np.empty(

R.<x,y,z> = QQ[]
I = R*(x+y*z,2*x^2+x*z+y+1)
F = I.gens()
#x0 = np.array([[-1,0],[0,1]]) # puiseux initial form solution
x0 = np.array([[-1,0,-3,0],[0,1,0,-3]]) # "" with one more term found
asdf = newtonStep(x0,F,jacobian(F,(y,z)),termsTarget=4)
print asdf

#------------------------------------------------------------------------------#
def alternateNewtonStep(ci,x0evalled,F,J,termsTarget,recompute=False):
    """
    Same as above, except we return x1 evaluated at the ci's plus
    at new points, coupled with the new list of points with our new ones
    appended.
    Computes next set of terms (coeffs of series).
    x0 is an nxk matrix, n=#variables, k=#terms, a guess for the root.
    F is a list of polynomials.
    J is the Jacobian of F.
    termsTarget is the number of terms to attempt to compute. Becomes k.
        Note: termsTarget shouldn't exceed #terms(x0), probably?
    recompute is whether we re-interpolate for the coeffs of x0 or just reuse
        them.
    """
    pass
