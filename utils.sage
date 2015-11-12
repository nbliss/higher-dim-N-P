"""
Computes the UCT based on Vector, which
is a vector or list of vectors.
"""
def uct(Vector):
    def ConvertMatrixToList(Matrix):
        Pts = []
        for Row in Matrix.transpose():
            Pts.append(list(Row))
        return Pts
    #Stack points vertically
    HNF, UCT = (matrix(matrix(Vector),ZZ).transpose()).echelon_form(include_zero_rows = True, transformation = True)
    #We need to check to make sure we're getting the right
    UCT = ((UCT)^-1)
    return UCT.transpose()

"""
Given a matrix mat (from a UCT) and a ring R,
constructs a dictionary giving the change-of-variables
associated with mat.
"""
def changeVariables(I,mat,R):
    subDict = {}
    for i in xrange(len(mat.columns())):
        col = mat.columns()[i]
        toAdd = 1
        for j in xrange(len(col)):
            toAdd *= R.gens()[j]^col[j]
        subDict[R.gens()[i]] = toAdd
    #return I.subs(in_dict = subDict),subDict
    return subDict

def initialForm(p,v):
    """
    Compute the initial form of p wrt v
    """
    v = vector(v)
    d = p.dict()
    exps = d.keys()
    smallestExps = [exps[0]]
    minDot = vector(exps[0])*v
    for exponent in exps[1:]:
        m = vector(exponent)*v
        if m<minDot:
            minDot = m
            smallestExps = [exponent]
        elif m==minDot:smallestExps.append(exponent)
    return p.parent()({e:d[e] for e in smallestExps})
