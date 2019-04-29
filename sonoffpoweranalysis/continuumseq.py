import numpy, pandas
import scipy.optimize

eps0 = 3
eps1 = eps0/2

def quadraline(xs, a, xc, ya):
    return a*numpy.minimum(0, xs-xc)**2 + ya

class continuumseq:
    def __init__(self, vs, i0, i1, cseqs):
        self.vs = vs
        self.i0, self.i1 = i0, i1
        if cseqs:
            self.prevcseq = cseqs[-1]
            cseqs[-1].nextcseq = self
        else:
            self.prevcseq = None
        self.nextcseq = None
        self.initfixedvals()
        self.initfixedqvals(eps1)
        
    def initfixedvals(self):
        self.leng = self.i1 - self.i0
        self.cvs = self.vs[self.i0:self.i1]
        self.mean = self.cvs.mean()
        self.sqsum = self.cvs.var()*self.leng  # sum((cseq.cvs - cseq.mean)**2)
        self.clo, self.chi = min(self.cvs), max(self.cvs)

    def initfixedqvals(self, eps1):
        self.a = 0
        if self.leng >= 5 and self.clo < self.chi:
            xs = numpy.arange(self.leng)
            k = scipy.optimize.curve_fit(quadraline, xs, self.cvs, p0=(0, self.leng, self.mean), bounds=([0, 0, self.clo], [1, self.leng, self.chi]))
            self.a, self.xc, self.ya = k[0]
            if (self.a*self.xc**2) < eps1:
                self.a = 0  # quadratic section so minor not worth it
            else:
                self.sqsumq = sum((self.cvs - quadraline(xs, self.a, self.xc, self.ya))**2)
        if self.a == 0:
            self.a, self.xc, self.ya = 0, 0, self.mean
            self.sqsumq = self.sqsum
        
    def quadraline(self, xs):
        return quadraline(xs, self.a, self.xc, self.ya)
        
    def getqcvs(self):
        xs = numpy.arange(self.leng)
        return quadraline(xs, self.a, self.xc, self.ya)
    
    def getflatcvs(self):
        return self.cvs[int(self.xc):]
    
    
    
# find the continuum segments
def makecontinuumsegments(vs):
    cseqs = [ ]
    i0 = 0
    for i in range(1, len(vs)):
        if abs(vs[i] - vs[i-1]) >= eps0:
            cseqs.append(continuumseq(vs, i0, i, cseqs))
            i0 = i
    cseqs.append(continuumseq(vs, i0, len(vs), cseqs))
    return cseqs


