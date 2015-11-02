def generalNewtonPuiseuxTest():
    load("generalNewtonPuiseux.sage")
    print 'hi'*30
    R.<x,y,z> = PolynomialRing(QQ,3)
    p = x^2 + y^2 + z^2 + 4*x
    q = x^2 + y^2 + 2*x
    I = (p,q)*R
    SOLUTION = newtonPuiseux(I)

def pSeriesTupleTest():
    load("pSeriesTuple.sage")
    ps = pSeriesTuple()
    print ps
    ps.addTerm([2,2,2],[1,1,1])
    print ps
    ps.addTerm([2,1,1],[3,3,3])
    print ps

def trop_intersection_wrapperTest():
    load("trop_intersection_wrapper.sage")
    R.<x,y,z> = PolynomialRing(CC,3)
    p = x^2 + y^2 + z^2 + 4*x
    q = x^2 + y^2 + 2*x
    I = (p,q)*R
    print getInitialForms(I)

testDict = {'generalNewtonPuiseux':generalNewtonPuiseuxTest,'pSeriesTuple':pSeriesTupleTest,'trop_intersection_wrapper':trop_intersection_wrapperTest}

import sys
tests = []
try: 
    choice = str(sys.argv[1])
    if '.sage' in choice:choice = choice.split('.')[0]
    if choice=='all':
        tests = testDict.values()
        tests.remove(generalNewtonPuiseuxTest)
    else:tests = [testDict[choice]]
except Exception:
    print "invalid entry"

for c in tests:
    c()
    print '='*30
