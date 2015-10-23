load("utils.sage")
"""
Outline:
At each (recursively done) step, program asks user for which
cone to continue with
"""


def getRationalCoeffs(I,clockout = 3,height_bound = 0):
    """
    Searches for a rational solution to the ideal defined by I.
    If height_bound is specified and greater than 0, uses it as
    the bound, otherwise increments the bound until points are found
    or clockout time is spent, whichever comes first.
    Clockout doesn't quite work--it will start a rational_points
    calculation as long as we haven't reached it, meaning it could
    start one just before the clockout and then take a while.
    WARNING: definitely might return an empty list even if variety
    is non-empty and contains rational points!!!
    """
    A = AffineSpace(I.ring()).subscheme(I)
    if height_bound > 0:
        points = A.rational_points(bound=height_bound)
        return filter(lambda pt:pt[0]==points[0][0],points)
    from time import time
    startTime = time()
    i = 1
    while time()-startTime<=clockout:
        points = filter(lambda pt:pt[0]!=0,A.rational_points(bound=i))
        if points==[]: i+=1
        else:
            return filter(lambda pt:pt[0]==points[0][0],points)


def getCoeffs(I):
    smallring = I.ring().remove_variable(R.0)
    v = (smallring*I.subs(x=1)).variety()
    toReturn = []
    for pt in v:
        toAdd = [1]
        toAdd+[pt[i] for i in R.gens()]
        toReturn.append(toAdd)
    return toReturn



EXPONENTS = []
COEFFS = []
RAMIND = 0

def maintainSeriesList(v,c):
    global RAMIND,EXPONENTS,COEFFS
    if EXPONENTS ==[]:
        EXPONENTS.append(v[1:])
        COEFFS.append(c)
        RAMIND = v[0]
    else:
        assert v[0]==1
        oldExp = EXPONENTS[-1]
        EXPONENTS.append([v[i+1]+oldExp[i] for i in xrange(len(v)-1)])
        COEFFS.append(c)

def npSubstitution(I,exps,coeffs):
    """
    Finds the next ideal using exps and coeffs,
    i.e. does the higher-dim version of substituting
    x^gamma(c+y) for y as we would do in 2 dimensions
    """
    subDict = {}
    R = I.ring()
    for i in xrange(1,len(exps)):
        thisGen = R.gens()[i]
        subDict[thisGen] = (R.0)^exps[i]*(coeffs[i]+thisGen)
    subDict[R.0] = coeffs[0]*(R.0)^exps[0]
    subbedIdeal = I.subs(subDict)
    return subbedIdeal


# takes an ideal I
def newtonPuiseux(I):
    if input("type y if done: ")=='y':
        var('t')
        ystuff = sum(COEFFS[i][1]*t^EXPONENTS[i][0] for i in xrange(len(COEFFS)))
        zstuff = sum(COEFFS[i][2]*t^EXPONENTS[i][1] for i in xrange(len(COEFFS)))
        answer = [COEFFS[0][0]*t^RAMIND]
        answer += []
        return answer
    R = I.ring()
    S = LaurentPolynomialRing(R.base_ring(),R.variable_names())
    F = I.groebner_fan().tropical_intersection()
    inForms = F.initial_form_systems()
    i = 0
    for form in inForms:
        print 'i='+str(i); i+=1
        rays = 0-matrix(form.rays())
        subDict = changeVariables(form.initial_forms()*R,uct(rays),S)
        print "Cone: "
        print rays
        print "Initial form: ",[factor(f) for f in form.initial_forms()]
        print "With x=1: ",[f.subs(x=1) for f in form.initial_forms()]
        """
        print "Substitution from UCT: ",subDict
        sdf = [factor(S(f).subs(in_dict = subDict)) for f in form.initial_forms()]
        print "Post-substitution: ",sdf
        print "Without units: ",[expand(f/f.unit()) for f in sdf]
        """
        print
    print '-'*35
    toExpand = int(input("Choose a cone using \'i\'"))
    form = inForms[toExpand]
    v = [0-a for a in form.rays()]
    c = getCoeffs(form.initial_forms())*R
    i = 0
    for pt in c:
        print 'i=',i
        i+=1
        print pt
    c = c[int(input("Choose coeff: "))]
    maintainSeriesList(v,c)
    
    newtonPuiseux(npSubstitution(I,v,c))
