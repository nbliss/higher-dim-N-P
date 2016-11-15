load("puiseuxNewtonMethod.sage")
load("../tropicalBasisStuff.sage")
load("../tropVarComputations.sage")
load("../generalNewtonPuiseux.sage")
def getStartVector(ringVar,I,V):
    initSyst = initialGeneratorIdeal(I,V)
    c = getRationalCoeffs(initSyst)
    if c==[]:c = getCoeffs(initSyst)
    if len(c)>1:
        i = 0
        for pt in c:
            print 'i='+str(i)+':'
            i+=1
            print pt
        c = c[getInput("Choose coeff by giving i--> ",int)]    
    else:c=c[0]
    return vector([coeff*ringVar^v for coeff,v in zip(c,V)])
    
def rayHomotopy(I,V,varIndex = 0):
    """
    Runs a homotopy on I from a general system to the ideal I.
    Tracks the series starting with ray V.
    """
    R = I.ring()
    startSyst = randomCoeffs(I)
    targetSyst = makeIdealOverPS(I,varIndex=varIndex)[0]
    R_PS = targetSyst.ring()
    x0 = getStartVector(R_PS.gen(0),startSyst,V)
    startSyst = makeIdealOverPS(startSyst,varIndex=varIndex)[0]
    Rt = PolynomialRing(R_PS,'t')
    t = Rt.0
    startSyst,targetSyst = map(lambda a:Rt*a,[startSyst,targetSyst])
    homotopySyst = Rt*[(1-t)*p1 + t*p2 for p1,p2 in zip(startSyst.gens(),targetSyst.gens())]
    tval = 0
    F = lambda tempTVal:R_PS*map(lambda a:R_PS(a.subs(tempTVal)),homotopySyst.gens())
    print '1--------------------'*4
    nextStep = newtonMethod(x0,F(tval))
    print '2--------------------'*4
    intermediateSols = [nextStep]
    tval = tval + 0.05
    intermediateSols.append(newtonMethod(intermediateSols[-1],F(tval)))
    for i in xrange(19):
        tval = tval + 0.05
        nextGuess = 2*intermediateSols[-1] - intermediateSols[-2]
        intermediateSols.append(newtonMethod(nextGuess,F(tval)))
    return intermediateSols[-1]

example = 1
if example==1:
    """
    target system ray: [1,1,1]
    start system ray: [1,3,3] (generic coefficient system)
    """
    R.<x,y,z> = QQ[]
    targetSyst = R*[y^2 - x^3 + y + z, z^2 + x^3 + y + z]
    asdf = rayHomotopy(targetSyst,[1,3,3])
    print asdf
