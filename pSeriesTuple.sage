class pSeriesTuple(object):
    def __init__(self,expander_index=0):
        self.coeffs = []
        self.exps = []
        self.expander_index = expander_index
        self.numVars = 0

    def addTerm(self,coeffs,exps):
        self.coeffs.append([c for c in coeffs])
        if type(exps)!=list:
            exps = list(exps[0])
        if self.exps==[]:
            self.exps.append([e for e in exps])
            self.numVars = len(coeffs)
            return
        toAdd = [e for e in exps]
        for i in xrange(self.numVars):
            if i==self.expander_index: continue
            else: toAdd[i] += self.exps[-1][i]
        self.exps.append(toAdd)

    def seriesTuple(self):
        if self.exps==[]:return []
        myRing.<t> = self.coeffs[-1][0].base_ring()[]
        toReturn = [0]*self.numVars
        #   a = a1  a2  a3  a4  (coeffs)
        #   g = g1  g2  g3  g4  (exps)
        #
        #    ((a4^g3 * a3)^g2 * a2)^g1 * a1
        for variableNum in xrange(self.numVars):
            if variableNum==self.expander_index:
                expander_coeff=1
                expander_exp = 1
                for i in xrange(len(self.coeffs)):
                    expander_exp *= self.exps[i][self.expander_index]
                    i = len(self.coeffs)-i-1
                    expander_coeff*=self.coeffs[i][self.expander_index]
                    if i-1>=0:expander_coeff**=self.exps[i-1][self.expander_index]
                toReturn[self.expander_index] = expander_coeff * t ** expander_exp
                continue
            for termNum in xrange(len(self.coeffs)):
                thisCoeff = self.coeffs[termNum][variableNum]
                thisExp = self.exps[termNum][variableNum]
                toReturn[variableNum] += thisCoeff * t ** thisExp
        return toReturn
    def __eq__(self,other):
        if type(other)!=pSeriesTuple or self.numVars != other.numVars:
            return False
        for i in xrange(self.numVars):
            if self.seriesTuple()[i]!=other.seriesTuple()[i]:return False
        return True
    def __repr__(self):
        tup = self.seriesTuple()
        baseRing = max([coeff.base_ring() for coeff in self.coeffs[-1]])
        myRing = PowerSeriesRing(baseRing,'t')
        return "Tuple "+str([myRing(elt) for elt in tup])+" of Puiseux series"
    def __str__(self):return self.__repr__()
