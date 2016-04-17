"""
Functions for computing tropical stuff.
Mostly consists of algorithms taken from Computing Tropical
Varieties, or code that calls gfan.
"""
load("tropicalBasisStuff.sage")
load("badPrevars.sage")

def mCutDownCone(I,vects):
    """
    Interfaces to macaulay2 to compute a witness
    that rules out vects from the tropical prevariety,
    if one exists.
    """
    R = I.ring()
    macaulay2('QQ['+','.join(R.variable_names())+']')
    macaulay2('load "mCutDownCones.m2"')
    macaulay2('I = ideal('+str(I.gens())[1:-1]+')')
    vectString = str(vects).replace('[','{').replace(']','}')
    returned = macaulay2("cutDownCone(I,"+vectString+")")
    return returned.to_sage()

def mTropBasis(I,F=None,dim=None):
    """
    dim is an int. If given, we only cut down
    cones of that dimension.
    """
    if F==None:F = getPrevar(I)
    polys = []
    def cutDimN(i):
        for cone in F.cones()[i]:
            rays = [F.rays()[j] for j in cone]
            if rays==[]:return
            rays = [[-a for a in v] for v in rays]
            p = mCutDownCone(I,rays)
            if p!=0:
                polys.append(p)
    if dim!=None and dim<=F.dim():cutDimN(dim)
    else:
        for i in xrange(1,F.dim()+1): cutDimN(i)
    return I if polys==[] else I+polys

def mTropVar(I,F=None):return getPrevar(mTropBasis(I,F))

def monInIdeal(I):
    """
    Returns a monomial, if one exists, in I,
    else returns 0.
    """
    R = I.ring()
    #m = x1*x2*...*xn
    m = reduce(lambda a,b:a*b,R.gens())
    varProd = m
    if I.saturation(m)[0]!=R*1:return False
    while m not in I: m *= varProd
    return m

def cutDownCone(I,vects):
    """
    For a cone (list of lists) "vects" that SHOULD
    be in the prevariety, returns 0 if that cone is
    in the variety or 'p' such that I+<p> doesn't contain
    the cone "vects" in its tropical intersection.
    """
    innerRay,otherguy = getRay(vects,I)
    #print innerRay
    I = I.homogenize()
    R = I.ring()
    #weightedR = PolynomialRing(R.base_ring(), R.variable_names(), order=TermOrder('negwdegrevlex',innerRay))
    weightedR = PolynomialRing(R.base_ring(), R.variable_names(), order=TermOrder('wdeglex',innerRay))
    weightedI = weightedR*I
    xm = monInIdeal(initialIdeal(I,innerRay))
    #print innerRay
    #print initialIdeal(I,innerRay)
    if xm==0:return 0
    f = xm - weightedI.reduce(xm)
    f = f.subs({R.gens()[-1]:1})
    #assert f in I
    #return f,otherguy[:2]
    return f#,xm


def getRay(vects,I):
    if type(vects[0])==list:
        if len(vects)>1:
            innerRay = list(sum([randint(1,1000)*vector(v) for v in vects]))
        else:innerRay = vects[0]
    else:innerRay = vects
    #print initialIdeal(I,innerRay)
    oldy = [i for i in innerRay]
    #innerRay = [1,3,1,1]
    #print innerRay,'-'*20
    maxy = max(innerRay)+1
    innerRay = [maxy-i for i in innerRay]+[maxy]
    #print initialIdeal(I.homogenize(),innerRay).subs({4:1})
    #return [841, 1, 841, 841, 975], oldy
    #return [2,1,2,2,3],oldy
    return innerRay,oldy

    newVects = []
    for v in vects:
        maxy = max(v)+1
        newVects.append([maxy-i for i in v]+[maxy])
    vects = newVects
    #print vects
    if type(vects[0])==list:
        if len(vects)>1:
            innerRay = list(sum([randint(1,1000)*vector(v) for v in vects]))
        else:innerRay = vects[0]
    else:innerRay = vects
    #maxy = max(innerRay)+1
    #innerRay = [maxy-i for i in innerRay]+[maxy]
    print innerRay
    return innerRay


def OTHERcutDownCone(I,vects):
    """
    For a cone (list of lists) "vects" that SHOULD
    be in the prevariety, returns 0 if that cone is
    in the variety or 'p' such that I+<p> doesn't contain
    the cone "vects" in its tropical intersection.
    """
    if type(vects[0])==list:
        if len(vects)>1:
            innerRay = list(sum([randint(1,1000)*vector(v) for v in vects]))
        else:innerRay = vects[0]
    else:innerRay = vects
    R = I.ring()
    weightedR = PolynomialRing(R.base_ring(), R.variable_names(), order=TermOrder('negwdegrevlex',innerRay))
    weightedI = weightedR*I
    xm = monInIdeal(initialIdeal(weightedI,innerRay))
    if xm==0:return 0
    f = xm - weightedI.reduce(xm)
    #assert f in I
    return f

def initialIdeal(I,v,useGfan = True):
    """
    Computes the initial ideal of I wrt v.
    Fails if v has any negative entries.
    """
    if not useGfan:
        # Good chance this is wrong!!!!!!
        print 'This initial ideal calculation might be wrong!!!!!!'
        R = PolynomialRing(I.base_ring(), I.ring().variable_names(), order=TermOrder('negwdegrevlex',v))
        gb=R*(R*I).groebner_basis()
        #return gb
        return I.ring()*initialGeneratorIdeal(gb,v)

    # -----Gfan try------ #
    R = I.ring()
    if not I.is_homogeneous():
        newVar = 'w'
        hI = (PolynomialRing(R.base_ring(),R.variable_names())*I).homogenize(newVar)
        #hR = PolynomialRing(R.base_ring(),[newVar]+list(R.variable_names()))
        #hI = I.homogenize(newVar)
    else:
        hI = I
    hR = hI.ring()
    filestring = "Q["+reduce(lambda a,b:a+','+b,hR.variable_names())+']\n'
    filestring+='{' + str(hI.gens())[1:-1] + '}\n'
    if type(v) not in [list,tuple]:v = list(v[0])
    if not I.is_homogeneous():
        filestring+='\n'+str([-vi for vi in v]+[0])
    else:filestring+='\n'+str([-vi for vi in v])
    inputfile = '/tmp/myfile'
    outputfile = '/tmp/gfanOut'
    f = open(inputfile,'w')
    f.write(filestring)
    f.close()
    os.system('/Applications/gfan_files/gfan_initialforms --ideal < '+ inputfile+' > '+outputfile)
    g = open('/tmp/gfanOut')
    toReturn = g.read()
    g.close()
    os.remove(inputfile)
    os.remove(outputfile)
    if not I.is_homogeneous():
        toReturn = toReturn.replace('w','1')
    toReturn = toReturn.split('{')[1][:-2].split(',')
    #toReturn = hR*[hR(poly[1:]) for poly in toReturn]
    #toReturn = R*toReturn.subs(in_dict={hR.gens()[-1]:1})
    toReturn = R*[R(poly[1:]) for poly in toReturn]
    return toReturn

    # -----Second try------ #
    # Trying to do negative powers, but didn't really work
    """
    if reduce(lambda a,b:a and b>0,v,True) or reduce(lambda a,b:a or b<0,v,False):
        R = PolynomialRing(I.base_ring(), I.ring().variable_names(), order=TermOrder('negwdegrevlex',v))
        return I.ring()*initialGeneratorIdeal(R*(R*I).groebner_basis(),v)
    else:
        oldRing = I.ring()
        zeroVars = []
        posVars = []
        for i in xrange(len(v)):
            thisVar = oldRing.variable_names()[i]
            if v[i]==0:zeroVars.append(thisVar)
            else:posVars.append(thisVar)
        baseRing = PolynomialRing(I.base_ring(),zeroVars)
        #baseRing.inject_variables()
        R = PolynomialRing(baseRing,posVars, order=TermOrder('negwdegrevlex',filter(lambda a:a>0,v)))
        #R.inject_variables()
        gb = (R*[R(a._repr_()) for a in I.gens()]).groebner_basis('toy:buchberger')
        return initialGeneratorIdeal(R*gb,v)
    """

def getPrevar(I):
    return I.groebner_fan().tropical_intersection()

def initialForm(p,v):
    """
    Computes the initial form of p wrt v
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

def initialGeneratorIdeal(I,v):
    return I.ring()*map(lambda a:initialForm(a,v),I.gens())
