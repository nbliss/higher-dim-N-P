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
    N = R.ngens()-1
    print "x0:\n",x0
    k = x0.shape[1] # current number of known terms
    ci = np.random.rand(termsTarget,1)/1000
    # fixed seed for predictability:
    #ci = np.array([[i*0.1] for i in xrange(termsTarget)])
    #ci = np.array([[i] for i in xrange(2,termsTarget+2)])
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
    LHS = []
    RHS = []
    # Evaluating J and F, and multiply on LHS by J.
    # These were horizontal slices in my notes.
    for i in xrange(termsTarget):
        subDict = {xi: x0ci for xi,x0ci in zip(R.gens()[1:],x0evalled[i])}
        subDict[R.0] = ci[i][0] # for first variable, sub just the value
        Ji = np.array(map(lambda a:a.subs(subDict),J.list()))
        Ji.shape = (N,N)
        Fi = -np.array(map(lambda f:[f.subs(subDict)],F)) #col vector of F@ci
        LHS.append(Ji)
        RHS.append(Fi)   # List of: #vars x 1 matrices
    RHS = np.concatenate(RHS) #vars x 1 x #ci
    
    leftC = np.insert(powers,range(1,termsTarget+1),0,axis=0)
    rightC = np.insert(powers,range(termsTarget),0,axis=0)
    C = np.concatenate([leftC,rightC],axis=1)

    zeros = np.zeros(LHS[0].shape)
    for i in xrange(len(LHS)):
        asdf = [zeros]*termsTarget
        asdf[i] = LHS[i]
        LHS[i] = np.concatenate(asdf,axis=1)
    Jmat = np.concatenate(LHS)
    LHS = np.dot(Jmat,C)
    sol = np.linalg.solve(LHS,RHS)
    sol.shape=[2,-1]
    x0 = np.insert(x0,[k]*(termsTarget-k),0,axis=1)
    return x0 + sol
    return LHS,RHS,sol


#------------------------------------------------------------------------------#

R.<x,y,z> = QQ[]
I = R*(x+y*z,2*x^2+x*z+y+1)
F = I.gens()
#x0 = np.array([[-1,0],[0,1]]) # puiseux initial form solution
x0 = np.array([[-1,0,-3,0],[0,1,0,-3]]) # "" with one more term found
asdf = newtonStep(x0,F,jacobian(F,(y,z)),termsTarget=12)
print "Solution:"
print asdf

"""
print "Solution:"
for i in asdf:
    print '[',
    for j in i[0]:print j,' ',
    print ']'
"""
