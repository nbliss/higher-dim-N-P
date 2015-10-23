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
