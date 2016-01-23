load("utils.sage")
load("trop_intersection_wrapper.sage")
load("pSeriesTuple.sage")
"""
Outline:
At each (recursive) step, program asks user for which
cone to continue with.
"""

#-----------------------------------------------------------------------#
def reducePoly(p):
    """
    Factors out any extra x_i's.
    UNLESS p is a monomial
    """
    if p==0:return p
    exps = p.exponents()
    if len(exps)==1:return p
    R = p.parent()
    nVars = R.ngens()
    toSubtract = [0]*nVars
    for i in xrange(nVars):
        toSubtract[i] = min(map(lambda a:a[i],exps))
    d = p.dict()
    newDict = {}
    for exponent in p.exponents():
        newDict[tuple([int(exponent[i]-toSubtract[i]) for i in xrange(nVars)])] = d[exponent]
    return R(newDict)

def reduceIdeal(I):
    """
    Returns ideal with extra x_i's factored out of generators
    """
    return I.ring()*map(lambda a:reducePoly(a),I.gens())

#-----------------------------------------------------------------------#
def getRationalCoeffs(I,clockout = 2,height_bound = 0):
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
    I = reduceIdeal(I)
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
    I = reduceIdeal(I)
    R = I.ring()
    smallRing = R.remove_var(R.gens()[0])
    smallIdeal = (smallRing*I.subs({R.0:1}))
    v = smallIdeal.variety()
    if v==[]:v = smallIdeal.variety(CC)
    if v==[]:raise Exception("Initial form system has no solutions!!")
    newRing = v[0].keys()[0].parent()
    toReturn = []
    for pt in v:
        toAdd = [1]
        toAdd+=[pt[i] for i in newRing.gens()]
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
    f = form.initial_forms()
    print "initial form: ",f
    R = f[0].parent()
    #print "Initial form: ",[f.factor(proof=False) for f in form.initial_forms()]
    print "With x=1: ",[f.subs({R.0:1}) for f in form.initial_forms()]
    vol = form.mixedVolume()
    if vol!=None:print "mixed volume: ",vol
    """
    S = LaurentPolynomialRing(QQbar,R.variable_names())
    subDict = changeVariables(form.initial_forms()*R.change_ring(QQbar),uct(rays))
    print "Substitution from UCT: ",subDict
    #sdf = [factor(S(f).subs(in_dict = subDict),proof=False) for f in form.initial_forms()]
    sdf = [S(f).subs(in_dict = subDict) for f in form.initial_forms()]
    print "Post-substitution: ",sdf
    print "With x=1: ",[f.subs(x=1) for f in sdf]
    #print "Without units: ",[expand(f/f.unit()) for f in sdf]
    """
    uctSubbed = uctSub(R*form.initial_forms(),rays)
    print "Post-substitution from UCT: ",uctSubbed 
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
        print "Invalid entry \'%s\', expected %s. Type 'q' to quit." %(toReturn,myType)
        return getInput(s,myType)

#-----------------------------------------------------------------------#
# takes an ideal I
def performStep(I,SOLUTION):
    R = I.ring()
    print '-'*44
    print I
    """
    print "Linear portion:"
    for p in I.gens():
        linearExps = filter(lambda a:sum(a[1:])<2,p.dict().keys())
        print R({e:p.dict()[e] for e in linearExps})
    """
    """
    completeTest = I.subs(in_dict={i:0 for i in R.gens()[1:]})
    if completeTest==R*0:
        print 'Done!'
        return SOLUTION
    """
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
        print "Only one exponent possibility. "
        printConeStuff(form)
    else: 
        if len(inForms)==0:
            print "No satisfactory rays! Exiting..."
            return SOLUTION
            """
            print "No satisfactory rays! Printing all..."
            inForms = oldInForms
            """
        i = 0
        for form in inForms:
            print 'i='+str(i)+':'; i+=1
            printConeStuff(form)
        toExpand = getInput("Choose a cone by giving \'i\'--> ",int)
        form = inForms[toExpand]

    if form.rays().nrows()>1:
        print "You chose the higher dim cone"
        print form.rays()
        form = form.changeRays(getInput("Enter a ray inside it: ",eval))
        print form.rays()
    rational = 'n'
    #if SOLUTION.seriesTuple()==[]:
    rational = getInput("Try for rational coeffs? (y/anything else) ",str)
    c = []
    if rational=='y':
        #heightBound = 0
        #heightBound = getInput("Set height bound--> ",int)
        c = getRationalCoeffs(form.initial_forms()*R)
    if c==[]:
        if rational=='y':print "Sorry, no rational points found."
        c = getCoeffs(form.initial_forms()*R)
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
    if getInput("Type y if done: ",str)=='y':return SOLUTION
    #################################################################
    ###   Need to check if one of the series is finished!!!!!!!   ###
    #################################################################
    return performStep(npSubstitution(I,v,c),SOLUTION)


#-----------------------------------------------------------------------#
def newtonPuiseux(I):
    print 'Type \'q\' at any prompt to quit'
    SOLUTION = pSeriesTuple()
    return performStep(I,SOLUTION)

#-----------------------------------------------------------------------#
