class pSeriesTuple(object):
    def __init__(self,first_coeffs,first_exps,expander_index=0):
        assert len(first_coeffs)==len(first_exps)
        self.coeffs = [first_coeffs]
        self.exps = [first_exps]
        self.numVars = len(first_exps)
        self.expander_index = expander_index

    def addTerm(self,coeffs,exps):
        self.coeffs.append(coeffs)
        toAdd = exps
        for i in xrange(self.numVars):
            if i==self.expander_index: continue
            else: toAdd[i] += self.exps[-1][i]
        self.exps.append(toAdd)

    def seriesTuple(self,parameter = 't'):
        if type(parameter)==str:parameter = var(parameter)
        toReturn = [0]*self.numVars
        #   a = a1  a2  a3  a4
        #   g = g1  g2  g3  g4
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
                toReturn[self.expander_index] = expander_coeff * parameter ** expander_exp
                continue
            for termNum in xrange(len(self.coeffs)):
                thisCoeff = self.coeffs[termNum][variableNum]
                thisExp = self.exps[termNum][variableNum]
                toReturn[variableNum] += thisCoeff * parameter ** thisExp
        return toReturn

if __name__=="__main__":
    ps = pSeriesTuple([2,2,2],[1,1,1])
    print ps.seriesTuple()
    #ps.addTerm([1,1,1],[1,1,1])
    ps.addTerm([2,1,1],[3,3,3])
    print ps.seriesTuple()
