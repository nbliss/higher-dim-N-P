load("UCTutils.sage")
load("trop_intersection_wrapper.sage")
load("pSeriesTuple.sage")
"""
Outline:
At each (recursive) step, program asks user for which
cone to continue with.
"""


#-----------------------------------------------------------------------#
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
    points = []
    while time()-startTime<=clockout:
        points = filter(lambda pt:pt[0]!=0,A.rational_points(bound=i))
        if points==[]: i+=1
        else:
            return filter(lambda pt:pt[0]==points[0][0],points)
    return points


#-----------------------------------------------------------------------#
def getCoeffs(I):
    R = I.ring()
    smallRing = R.remove_var(R.gens()[0])
    smallIdeal = (smallRing*I.subs(x=1))
    v = smallIdeal.variety()
    if v==[]:v=smallIdeal.variety(ring=ComplexField(15))
    toReturn = []
    for pt in v:
        toAdd = [1]
        toAdd+=[pt[i] for i in smallRing.gens()]
        toReturn.append(toAdd)
    return toReturn


#-----------------------------------------------------------------------#
def npSubstitution(I,exps,coeffs):
    """
    Finds the next ideal using exps and coeffs,
    i.e. does the higher-dim version of substituting
    x^gamma(c+y) for y as we would do in 2 dimensions
    """
    if type(exps)!=list:exps = exps.list()
    subDict = {}
    R = I.ring()
    for i in xrange(1,len(coeffs)):
        thisGen = R.gens()[i]
        subDict[thisGen] = (R.0)^exps[i]*(coeffs[i]+thisGen)
    subDict[R.0] = coeffs[0]*(R.0)^exps[0]
    subbedIdeal = I.subs(subDict)
    return subbedIdeal

#-----------------------------------------------------------------------#
def printConeStuff(form):
    rays = matrix(form.rays())
    print "Cone: "
    print rays
    f = form.initial_forms()[0]
    #print "Initial form: ",[f.factor(proof=False) for f in form.initial_forms()]
    print "With x=1: ",[f.subs(x=1) for f in form.initial_forms()]
    """
    S = LaurentPolynomialRing(R.base_ring(),R.variable_names())
    subDict = changeVariables(form.initial_forms()*R,uct(rays),S)
    print "Substitution from UCT: ",subDict
    sdf = [factor(S(f).subs(in_dict = subDict)) for f in form.initial_forms()]
    print "Post-substitution: ",sdf
    print "Without units: ",[expand(f/f.unit()) for f in sdf]
    """
    print

#-----------------------------------------------------------------------#
def getInput(s,myType):
    """
    Gets input from user with prompt s, and coerces
    it to type myType, or repeates if failed.
    """
    toReturn = raw_input(s)
    if toReturn in ['q','Q']:
        import sys
        sys.exit()
    try:return myType(toReturn)
    except Exception:
        print "Invalid entry \'%s\'. Expected %s" %(toReturn,myType)
        return getInput(s,myType)

#-----------------------------------------------------------------------#
# takes an ideal I
def performStep(I,SOLUTION):
    R = I.ring()
    inForms = getInitialForms(I)
    oldInForms = [f for f in inForms]
    if SOLUTION.seriesTuple()==[]: #only want positive x exps for the first term
        def xPos(form):
            for ray in list(form.rays()):
                if ray[0] <= 0:return False
            return True
        inForms = filter(xPos,inForms)
    else: # want all exps >0 for subsequent terms
        def allPos(form):
            for ray in list(form.rays()):
                for element in ray:
                    if element <= 0:return False
            return True
        inForms = filter(allPos,inForms)
    if len(inForms)==1:
        form = inForms[0]
        print "Only one exponent possibility: ",form.rays()
    else: 
        if len(inForms)==0:
            print "No satisfactory rays! Printing all..."
            inForms = oldInForms
        i = 0
        for form in inForms:
            print 'i='+str(i)+':'; i+=1
            printConeStuff(form)
        print '-'*35
        toExpand = getInput("Choose a cone by giving \'i\'--> ",int)
        form = inForms[toExpand]
    rational = getInput("Try for rational coeffs? (y/anything else) ",str)
    c = []
    if rational=='y':
        heightBound = 0
        #heightBound = getInput("Set height bound--> ",int)
        c = getRationalCoeffs(form.initial_forms()*R,heightBound)
    if c==[]: c = getCoeffs(form.initial_forms()*R)
    v = form.rays()
    if len(c)==1:
        c = c[0]
        print "Only one coefficient possibility: ",c
    else:
        i = 0
        for pt in c:
            print 'i='+str(i)+':'
            i+=1
            print pt
        c = c[getInput("Choose coeff by giving i--> ",int)]
    SOLUTION.addTerm(c,v)
    print SOLUTION
    if getInput("type y if done: ",str)=='y':return SOLUTION
    return performStep(npSubstitution(I,v,c),SOLUTION)


#-----------------------------------------------------------------------#
def newtonPuiseux(I):
    print 'Type \'q\' at any prompt to quit'
    SOLUTION = pSeriesTuple()
    return performStep(I,SOLUTION)

#-----------------------------------------------------------------------#
