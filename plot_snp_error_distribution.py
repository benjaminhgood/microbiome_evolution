import numpy
import pylab
from scipy.stats import poisson
from scipy.stats import binom

Ds = numpy.array([20.0*0.3, 10,20,40])
dfs = numpy.linspace(0,1,11)

for D in Ds:

    perrs = []
    for df in dfs:
        
        perrs.append(binom.cdf((1-df)/2*D, D, 0.5)**2)
        
    pylab.semilogy(dfs, perrs,'.-',label=('D=%d' % D))

pylab.legend(loc='lower left',frameon=False)  

pylab.figure()

dfs = numpy.array([0.1,0.6,0.7,0.8])
Ds = numpy.arange(long(20*0.3),20)

for df in dfs:

    perrs = []
    for D in Ds:
        
        perrs.append(binom.cdf((1-df)/2*D, D, 0.5)**2)
        
    pylab.semilogy(Ds, perrs,'.-',label=('df=%g' % df))

pylab.legend(loc='lower left',frameon=False)  
  
pylab.show()
