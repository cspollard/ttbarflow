import numpy

MCCOEFFS = [ 4.27275009e-07 , -1.23074969e-02 , -7.58829585e-01 ]
DATACOEFFS = [ 5.59315967e-07 , -1.02652000e-02 , -1.59489289e+00 ]

def polyexp(coeffs, xs):
  return numpy.exp(numpy.polyval(coeffs, xs))

def MCfunc(xs):
  return polyexp(MCCOEFFS, xs)


def datafunc(xs):
  return polyexp(DATACOEFFS, xs)
