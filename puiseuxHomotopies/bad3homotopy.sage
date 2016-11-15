def F(syst,tval):
    """
    Return the system at a given t-point in time, over
    the correct ring (i.e. without t).
    Assumes ring is in variable t with coeffs over x,y,z, etc.
    INPUT: an ideal syst, a number tval
    OUTPUT: an ideal with tval substituted in, over the smaller
    ring.
    """
    ringWithoutT = (syst.0).coefficients()[0].parent()
    return ringWithoutT * map(lambda a:a.subs(tval),syst.gens())

"""
target system ray: [1,1,1]
start system ray: [1,3,3] (generic coefficient system)
"""
load("../tropicalBasisStuff.sage")
load("puiseuxNewtonMethod.sage")
R.<x,y,z> = QQ[]
targetSyst = R*[y^2 - x^3 + y + z, z^2 + x^3 + y + z]
#startSyst = randomCoeffs(targetSyst)
# vvv  for replicable runs
startSyst = R*(5/78*x^3 - 79/200*y^2 + 237/62*y - 57/50*z, -3/47*x^3 + 73/82*z^2 - 102/131*y + 37/292*z)

PP = PowerSeriesRing(ComplexField(120),'x',default_prec=25)
#PP = PowerSeriesRing(RealField(),'x',default_prec=15,sparse=True)
PP = PowerSeriesRing(RealField(),'x',default_prec=15)
#PP = PowerSeriesRing(QQ,'x',default_prec=35)
PP.inject_variables()
R.<y,z> = PP[]
R.inject_variables()
Rt.<t> = R[] # ring with t
Rt.inject_variables()
targetSyst = Rt*targetSyst
startSyst = Rt*startSyst 

homotopySyst = Rt*[(1-t)*startPoly + t*targetPoly for startPoly,targetPoly in zip(startSyst.gens(),targetSyst.gens())]
x0 = vector([1,-1])
tval = 0
nextStep = newtonMethod(x0,F(homotopySyst,tval))
intermediateSols = [nextStep]

tval = tval + 0.05
intermediateSols.append(newtonMethod(intermediateSols[-1],F(homotopySyst,tval)))

for i in xrange(19):
    tval = tval + 0.05
    nextGuess = 2*intermediateSols[-1] - intermediateSols[-2]
    intermediateSols.append(newtonMethod(nextGuess,F(homotopySyst,tval)))
