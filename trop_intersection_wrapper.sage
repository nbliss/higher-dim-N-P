"""
Wrapper for gfan's tropical computations. Allows us to compute
tropical prevarieties for ideals whose base ring is not QQ.
****** NOTE: RETURNS NEGATIVES OF GFAN'S RAYS, SINCE THIS IS WHAT
WE WANT FOR OUR ALGORITHM ******
"""
class inFormWrapper(object):
    def __init__(self,forms,rayList):
        self.rayList = rayList
        self.forms = forms
    def rays(self):return self.rayList
    def initial_forms(self):return self.forms

def getInitialForms(I):
    if I.base_ring()==QQ:
        systs = I.groebner_fan().tropical_intersection().initial_form_systems()
        return [inFormWrapper(f.initial_forms(),0-matrix(f.rays()))\
                for f in I.groebner_fan().tropical_intersection().initial_form_systems()]
    R = I.ring()
    polys = [p.dict() for p in I.gens()]
    qq_ring = I.ring().change_ring(base_ring = QQ)
    rationalCoeffPolys = []
    for i in xrange(len(polys)):
        rationalCoeffPolys.append(qq_ring(dict.fromkeys(polys[i],i+1)))
    # polys = [0.5x - y,  x^2 + 0.7 y]
    # rCP   = [  1x -1y, 2x^2 +   2 y]
    qq_ideal = qq_ring * rationalCoeffPolys
    F = qq_ideal.groebner_fan().tropical_intersection()
    forms = F.initial_form_systems()
    origForms = []
    for form in forms:
        origForm = []
        for p in form.initial_forms():
            d = p.dict()
            origPolyDict = polys[d.values()[0]-1]
            origForm.append(R({k:origPolyDict[k] for k in d.keys()}))
        origForms.append(inFormWrapper(origForm,form.rays()))
    return origForms
if __name__=='__main__':
    R.<x,y,z> = PolynomialRing(CC,3)
    p = x^2 + y^2 + z^2 + 4*x
    q = x^2 + y^2 + 2*x
    I = (p,q)*R
    print getInitialForms(I)
