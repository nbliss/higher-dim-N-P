"""
Wrapper for gfan's tropical computations. Allows us to compute
tropical prevarieties for ideals whose base ring is not QQ.
****** NOTE: RETURNS NEGATIVES OF GFAN'S RAYS, SINCE THIS IS WHAT
WE WANT FOR OUR ALGORITHM ******
"""
class inFormWrapper(object):
    def __init__(self,forms,rayList,rationalVersion = []):
        self.rayList = rayList
        self.forms = forms
        if rationalVersion==[]:
            self.rationalVersion = forms
        else:self.rationalVersion = rationalVersion
    def rays(self):return self.rayList
    def initial_forms(self):return self.forms
    def mixedVolume(self):
        R = self.rationalVersion[0].parent()
        smallRing = R.remove_var(R.gens()[0])
        theIdeal = smallRing*(R*self.rationalVersion).subs(x=1)
        try:return theIdeal.groebner_fan().mixed_volume()
        except Exception,e:print "Mixed volume failed. Message: %s" % str(e)

def getInitialForms(I):
    """
    If the ideal is already over QQ, no need for anything fancy. We still
    return our own inform wrapper for consistency, but do no more than call
    gfan. Otherwise, we convert the polynomials to dicts, construct new dicts
    with integer coefficients such that the coefficients are 1..n ordered to
    correspond to the order of the original polys in the original list. Pass
    these to gfan, then use the coeffs of the initial forms to figure out
    which of the original polys they correspond to, then make some new dicts
    and restore the original coefficents, finally wrapping as our objects.
    """
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
        origForms.append(inFormWrapper(origForm,0-matrix(form.rays()),form.initial_forms()))
    return origForms

