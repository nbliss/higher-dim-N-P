def badPrevar(n):
    """
    Returns an ideal in n variables with
    a bad 1-d cone.
    """
    varList = ['x'+str(i) for i in xrange(1,n+1)]
    R = PolynomialRing(QQ,varList)
    varList = R.gens()
    varSum = sum(varList)
    polyList = [- varList[0] + sum(varList[1:])+ varList[0]^2]
    for i in xrange(1,n-1):
        polyList.append(varSum + varList[i]^2)
    return R*polyList

def noetherBadPrevar(n):
    varList = ['x'+str(i) for i in xrange(1,n+1)]
    R = PolynomialRing(QQ,varList)
    varList = R.gens()
    varSum = sum(varList)
    polyList = [- varList[0] + sum(varList[1:])+ varList[1]^2]
    for i in xrange(1,n-1):
        polyList.append(varSum + varList[i+1]^2)
    return R*polyList


def otherBadPrevar(n):
    varList = ['x'+str(i) for i in xrange(1,n+1)]
    R = PolynomialRing(QQ,varList)
    varList = R.gens()
    varSum = sum(varList)
    polyList = [- varList[0] + sum(varList[1:])+ varList[0]^2]
    for i in xrange(1,n-1):
        polyList.append(varSum + varList[i]^2 + varList[i+1]^3)
    return R*polyList

