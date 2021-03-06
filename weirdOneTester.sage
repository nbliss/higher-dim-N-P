R.<x,y,z> = QQ[]
p1 = 4*x*z - y*z - z^2 + 3*x
p2 = 3*z^3 - 2/5*x*y - y*z - z^2 - 2/7*x
height_bound = 500
def randy():
    d1 = {a:QQ.random_element(height_bound,height_bound) for a in p1.dict().keys()}
    d2 = {a:QQ.random_element(height_bound,height_bound) for a in p2.dict().keys()}
    return R*(R(d1),R(d2))
def findTropVar(J):
    F = J.groebner_fan().tropical_basis()
    J = R*F
    return J.groebner_fan().tropical_intersection()

I2 = R*(43/246*x*z - 66/71*y*z + 287/402*z^2 - 88/169*x, -94/69*z^3 - 76/215*x*y + 367/315*y*z + 130/383*z^2 - 21/34*x)

n = 40
l = map(lambda a: a/n,range(1,n))

totalBad = 0
import sys
for i in xrange(500):
    print '\r'+str(i),
    sys.stdout.flush()
    asdf = findTropVar(randy()).rays()
    if asdf!=r:totalBad+=1
print
print 'height bound: ',height_bound
print totalBad


"""
500 trials

denominator bound 100: 21
denominator bound 300: 7
denominator bound 500: 6
"""
