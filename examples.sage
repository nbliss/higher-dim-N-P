load("generalNewtonPuiseux.sage")

def viviani():
    R.<x,y,z> = PolynomialRing(QQ,3)
    p = x^2 + y^2 + z^2 + 4*x
    q = x^2 + y^2 + 2*x
    I = (p,q)*R
    sol = newtonPuiseux(I)
    print sol

def curvyCircle():return
"""
R.<x,y,z> = PolynomialRing(QQ,3)
p = z^2+x^2+y^2-1
q = y^2-x^2-z
I = (p,q)*R
sol = newtonPuiseux(I)
print sol
"""

curvyCircle()
