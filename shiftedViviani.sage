load("utils.sage")

S.<x,y,z> = LaurentPolynomialRing(QQ,3)
R.<x,y,z> = PolynomialRing(QQ,3)
p = x^2 + y^2 + z^2 + 4*x
q = x^2 + y^2 + 2*x
I = (p,q)*R

def tropStep(I,verbose=False):
    if verbose:print '-'*35
    F = I.groebner_fan().tropical_intersection()
    inForms = F.initial_form_systems()
    if verbose:print "Finished? ",I.subs(y=0,z=0)==0
    i = 0
    for form in inForms:
        if verbose:print 'i='+str(i); i+=1
        rays = 0-matrix(form.rays())
        subDict = changeVariables(form.initial_forms()*R,uct(rays),S)
        if verbose:
            print "Cone: "
            print rays
            print "Initial form: ",[factor(f) for f in form.initial_forms()]
            print "With x=1: ",[f.subs(x=1) for f in form.initial_forms()]
            print "Substitution from UCT: ",subDict
        sdf = [factor(S(f).subs(in_dict = subDict)) for f in form.initial_forms()]
        if verbose:
            print "Post-substitution: ",sdf
            print "Without units: ",[expand(f/f.unit()) for f in sdf]
            print
    if verbose:print '-'*35
    return inForms

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

smallring.<y,z> = QQ[]
def getCoeffs(I):
    v = (smallring*I.subs(x=1)).variety()
    toReturn = []
    for pt in v:
        toReturn.append([1,pt[y],pt[z]])
    return toReturn

tropStep(I)


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

v=[2,1,1]
c=[-2,-2,-2] #roots are (-2, +/-2, +/-2)
maintainSeriesList(v,c)
I2 = npSubstitution(I,v,c)

for i in xrange(3):
    inForms = tropStep(I2)
    v = [0-a for a in inForms[0].rays()[0]] # [1,2,2]
    c = getCoeffs(inForms[0].initial_forms()*R)[0]
    maintainSeriesList(v,c)
    I2 = npSubstitution(I2,v,c)

var('t')
ystuff = sum(COEFFS[i][1]*t^EXPONENTS[i][0] for i in xrange(len(COEFFS)))
zstuff = sum(COEFFS[i][2]*t^EXPONENTS[i][1] for i in xrange(len(COEFFS)))
answer = (COEFFS[0][0]*t^RAMIND,ystuff,zstuff)
print answer
