import numpy
import scipy
from matplotlib import figure
from sys import argv
from cpplot import cpplot

MAXFITBINS = 8

toppt = \
  { "binning" : [300.0, 350.0, 400.0, 450.0, 500.0, 550.0, 650.0, 750.0, 1500.0]
  , "densities" : [0.00760485, 0.00485966, 0.00276359, 0.00158396, 0.00111556, 0.0006978, 0.000167224, 2.28222e-05]
  }

ttpt = \
  { "binning" : [0.0, 35.0, 75.0, 125.0, 170.0, 225.0, 800.0]
  , "densities" : [0.0100368, 0.00789039, 0.00333946, 0.00162792, 0.000813572, 8.36853e-05]
  }

ttm = \
  { "binning" : [225.0, 345.0, 420.0, 500.0, 590.0, 685.0, 790.0, 910.0, 1040.0, 1175.0, 1320.0, 3000.0]
  , "densities" : [0.000483109, 0.00281805, 0.0028295, 0.00207343, 0.0012649, 0.000762089, 0.000422668, 0.00022836, 0.00011654, 6.44363e-05, 7.16463e-06]
  }


def fitpolyexp(densities, binning):
  binning = numpy.array(binning)

  def polyexp(xs, a, b, c):
    coeffs = [a, b, c]
    return numpy.exp(numpy.polyval(coeffs, xs))

  centers = (binning[1:] + binning[:-1]) / 2

  optcoeff = scipy.optimize.curve_fit(polyexp, centers, densities, p0=numpy.zeros((3,)))[0]

  print(optcoeff)

  def f(xs):
    return polyexp(xs, optcoeff[0], optcoeff[1], optcoeff[2])

  return f


ptbinning = toppt["binning"]

MCinputs = numpy.loadtxt(argv[1], delimiter=",", skiprows=1)

toppthist = numpy.histogram(MCinputs[:,0], bins=ptbinning, density=True)[0]

MCfit = fitpolyexp(toppthist[:MAXFITBINS], ptbinning[:MAXFITBINS+1])
datafit = fitpolyexp(toppt["densities"], ptbinning[:MAXFITBINS+1])


fig = \
  cpplot.comparehist \
  ( [ cpplot.zeroerr(toppt["densities"]) , cpplot.zeroerr(toppthist) ]
  , numpy.array(toppt["binning"])
  , [ "data" , "Pythia8" ]
  , "$p_T^t$ / GeV"
  , "$\\frac{1}{\sigma}\\ \\frac{d\sigma}{dp_T^t}$ * GeV"
  , ratio=True
  )

plt = fig.axes[0]

xs = numpy.mgrid[ptbinning[0]:ptbinning[-1]:1000j]
mcys = MCfit(xs)
datays = datafit(xs)
plt.plot(xs, mcys, color="red", label="Pythia 8 fit")
plt.plot(xs, datays, color="gray", label="data fit")

plt.legend()
plt.set_yscale("log")

plt = fig.axes[1]
plt.set_ylabel("Pythia8 / data")
plt.plot(xs, mcys/datays, color="red")

fig.savefig("dataMC.png")

