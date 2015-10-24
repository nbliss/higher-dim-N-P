class inForm(object):
    def __init__(self,original_ring,original_polyDicts,gfan_inForm):
        self.rays = gfan_inForm.rays()
        polys = [p.dict() for p in gfan_inForms.initial_forms()]
        
        self.initial_forms = [original_ring(p) for 

    def rays(self):return self.rays
    def initial_forms(self):return self.initial_forms

def getInitialForms(I):
    polys = [p.dict() for p in I.gens()]
    qq_ring = I.ring().change_ring(base_ring = QQ)
    unitCoeffPolys = []
    for i in xrange(len(polys)):
        unitCoeffPolys.append(qq_ring(dict.fromkeys(polys[i],i)))
    qq_ideal = qq_ring*unitCoeffPolys
    F = qq_ideal.groebner_fan().tropical_intersection()
    inForms = F.initial_form_systems()

