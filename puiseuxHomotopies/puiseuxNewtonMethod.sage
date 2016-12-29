"""
x1 = x0 - J^-1(x0) * F(x0)
or
J(x0)*(x1-x0) = - F(x0)
using:
X = A.solve_right(B)
A*X == B

[a b] * F0(x0)
[c d]   F1(x0)
"""
def printy(a,rounded=False):
    allOfThem = map(lambda ind:zip(ind.coefficients(),ind.exponents()),a)
    N = max(map(len,allOfThem))
    for i in allOfThem:i += [('?','?')] * (N - len(i))
    allOfThem = map(list, zip(*allOfThem))
    for i in xrange(N):
        for c,e in allOfThem[i]:
            if type(c)==str:
                print '   \t\t',
                continue
            if rounded: print str(n(c,digits=4))+' x^'+str(e)+'\t',
            else: print str(c)+' x^'+str(e)+'\t',
        print

def newtonStep(x0,J,F):
    subDict = {a:b for a,b in zip(F[0].parent().gens(),x0.list())}
    J = J.subs(subDict)
    F = vector(map(lambda a:a.subs(subDict),F))
    #return (x0-J.solve_right(F))
    #return (x0-J.transpose().solve_left(F))
    return x0 - J.inverse()*F

def newtonMethod(sol,F,prec=15):
    """
    sol should be a vertical matrix
    """
    if type(F)==sage.rings.polynomial.multi_polynomial_ideal.MPolynomialIdeal:
        F = F.gens()
    J = jacobian(F,F[0].parent().gens())
    xCurr = newtonStep(sol,J,F)
    sols = []
    for i in xrange(15):
        xNext = newtonStep(xCurr,J,F)
        sols.append(xCurr)
        xCurr = xNext
    return sols[-1]

def makeIdealOverPS(I,varIndex=0,default_prec=35):
    origRing = I.ring()
    ringVars = list(origRing.gens())
    expandVar = ringVars.pop(varIndex)
    print ringVars
    #PP = PowerSeriesRing(origRing.base_ring(),str(expandVar),default_prec=default_prec,#sparse=True)
    PP = LaurentSeriesRing(origRing.base_ring(),str(expandVar),default_prec=default_prec)
    #PP = FractionField(origRing.remove_var(*ringVars))
    R = PolynomialRing(PP,str(ringVars)[1:-1])
    #R = PolynomialRing(origRing.remove_var(*ringVars))
    return R*I,PP

example=1
if example==1:
    """
    Viviani
    Puiseux series sol should be:
    y = 2*t - t^3 - 1/4*t^5 - 1/8*t^7 - 5/64*t^9 - 7/128*t^11 - 21/512*t^13,
    z =   2 - t^2 - 1/4*t^4 - 1/8*t^6 - 5/64*t^8 - 7/128*t^10 - 21/512*t^12
    """
    R.<x,y,z> = CC[]
    I = R*(y^2 + z^2 - 4 + 4*x^4, y^2 - 4*x^2 + 4*x^4)
    I,PP = makeIdealOverPS(I,default_prec=15)
    x=PP.0
    #x0 = vector([1.9999*x,1.9999])
    x0 = vector([2*x,2])
    asdf = newtonMethod(x0,I)
elif example==2:
    """
    First term isn't enough - regularity index is 2
    Puiseux series sol should be:
    y = -t^2 +   t^3 - 3*t^4 + 21/2*t^5 - 40*t^6 + 1275/8*t^7 - 1323/2*t^8 + 45193/16*t^9 - 24651/2*t^10 + 7004019/128*t^11
    z = -t^2 + 2*t^3 - 6*t^4 +   21*t^5 - 79*t^6 + 1263/4*t^7 -   1311*t^8 +  44789/8*t^9 -   24432*t^10 +  6942219/64*t^11
    """
    #R.<x,y,z> = ComplexField()[]
    R.<x,y,z> = QQ[]
    #I = R*(z^3 + x^2 - 2*x*y + y^2, y^3 + x - 2*y + z)
    I = R*(x^4 + 2*x^2*y + z^3 + y^2, y^3 - x^2 - 2*y + z)
    I,coeffRing = makeIdealOverPS(I)
    tvar = coeffRing.0
    #x0 = vector([-tvar^2,-tvar^2])
    x0 = vector([-tvar^2 + tvar^3,-tvar^2 + 2*tvar^3])
    #x0 = vector([-1,-1])
    #F = I.gens()
    #J = jacobian(F,I.ring().gens())
    #asdf = newtonStep(x0,J,F)
    asdf = newtonMethod(x0,F,prec=35)
elif example==3:
    """
    Hidden ray system.
    """
    PP = PowerSeriesRing(QQ,'x',default_prec=35)
    PP = PowerSeriesRing(ComplexField(100),'x',default_prec=35)
    PP.inject_variables()
    R.<y,z> = PP[]
    R.inject_variables()
    F = [y^2 - x^3 + y + z, z^2 + x^3 + y + z]
    jac = jacobian(F,[y,z])
    x0 = [-x,x]
    #x0 = [-x,1.1*x]
    asdf = newtonMethod(x0,F)
if example!=0:printy(asdf[-1],True)

