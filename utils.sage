def uct(Vector):
    """
    Computes the UCT based on Vector, which
    is a vector or list of vectors.
    """
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

def changeVariables(R,mat):
    """
    Given a matrix mat (from a UCT), constructs
    a dictionary giving the change-of-variables
    associated with mat.
    DEPRECATED: use uctSub instead
    """
    subDict = {}
    for i in xrange(len(mat.columns())):
        col = mat.columns()[i]
        toAdd = 1
        for j in xrange(len(col)):
            toAdd *= R.gens()[j]^col[j]
        subDict[R.gens()[i]] = toAdd
    return subDict

def uctSub(I,v):
    """
    Given a vector (or set thereof) v, finds the corresponding UCT matrix
    and makes the appropriate substitution in the ideal. Does so using
    vector mult., not the built-in, to avoid Sage's subst. order issues.
    
    """
    R = I.ring()
    mat = uct(v)
    newPolys = []
    for p in I.gens():   
        exps = mat*(matrix(p.exponents(as_ETuples=False)).transpose())
        exps = positivize(exps)
        oldExps = p.exponents()
        d = p.dict()
        newPoly = {tuple(exps.column(i)):d[oldExps[i]] for i in xrange(exps.ncols())}
        newPolys.append(R(newPoly))
    return R*newPolys

def positivize(mat):
    """
    Adds enough 1-vectors to each row to make them positive.
    Subtracts enough so that one row component is 0.
    """
    oneMatrix = matrix.ones(ZZ,mat.nrows(),mat.ncols())
    for i in xrange(mat.nrows()):
        m = min(mat.row(i))
        if m!=0:oneMatrix.set_row_to_multiple_of_row(i,i,-m)
        else:oneMatrix.set_row_to_multiple_of_row(i,i,0)
    return mat+oneMatrix
 
def support(I):
    """
    Return a list of support sets
    """
    toReturn = []
    for p in I.gens():
        toReturn.append(map(lambda a:list(a),p.exponents(as_ETuples=False)))
    return toReturn

"""
load("tropicalBasisStuff.sage")
def idealio(v):
    asdf = randomCoeffs(I,height_bound=1000)
    return initialGeneratorIdeal(asdf,v) == initialIdeal(asdf,v)

"""
def randomInitials():
    import sys
    for i in xrange(8,10):
        for j in xrange(1,10):
            for k in xrange(1,10):
                count = 0
                for numb in xrange(100):
                    if not idealio((i,j,k)):count += 1
                print str((i,j,k)),str(count),'\t',
                sys.stdout.flush()
