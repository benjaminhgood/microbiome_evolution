import numpy
import pylab
from scipy.stats import poisson

Li = 1e03
Lr = 1e02
c = 0.5
Dm = 20

fold_changes = numpy.array([2,3,4,5,6,7,8,9,10])

Navg = c*Dm*Li/Lr*1.0
Nmins = Navg/fold_changes

print Navg
print Nmins

perrs = poisson.cdf(Nmins, Navg)

pylab.semilogy(fold_changes, perrs,'b.-')
pylab.show()