import numpy
from scipy import optimize
from scipy import interpolate
from matplotlib import figure
from sys import argv
from cpplot import cpplot

# the binedge after which to begin the exponential polynomial fit
BINSPLIT = 5

# the degree of the exponential polynomial fit
POLYEXPDEG = 3


toppt = \
  { "binning" : [0.0, 25.0, 50.0, 75.0, 105.0, 135.0, 165.0, 195.0, 230.0, 265.0, 300.0, 350.0, 400.0, 450.0, 500.0, 1000.0]
  , "densities" : [0.00131799, 0.00382031, 0.0056759, 0.00627889, 0.0056098, 0.0040236, 0.00273743, 0.00179477, 0.00110479, 0.000751999, 0.000398752, 0.000217909, 0.000107683, 5.1266e-05, 7.13627e-06]
  }

ttpt = \
  { "binning" : [0.0, 35.0, 75.0, 125.0, 170.0, 225.0, 800.0]
  , "densities" : [0.0100368, 0.00789039, 0.00333946, 0.00162792, 0.000813572, 8.36853e-05]
  }

ttm = \
  { "binning" : [225.0, 345.0, 420.0, 500.0, 590.0, 685.0, 790.0, 910.0, 1040.0, 1175.0, 1320.0, 3000.0]
  , "densities" : [0.000483109, 0.00281805, 0.0028295, 0.00207343, 0.0012649, 0.000762089, 0.000422668, 0.00022836, 0.00011654, 6.44363e-05, 7.16463e-06]
  }


ptbinning = toppt["binning"]
centers = (numpy.array(ptbinning[:-1]) + numpy.array(ptbinning[1:]) ) / 2.0

if __name__ == "__main__":
  MCinputs = numpy.loadtxt(argv[1], delimiter=",", skiprows=1)
  toppthist = numpy.histogram(MCinputs[:,0], bins=ptbinning, density=True)[0]
  print(toppthist)


else:
  toppthist = [0.001120432548921514, 0.0034034365319120276, 0.005144407458005426, 0.00610997800055159, 0.005606331736288828, 0.0041407708144926975, 0.002944751239473296, 0.0019148041352732431, 0.0012643039937107782, 0.0008272470636096503, 0.00048536683919878396, 0.0002657251159300379, 0.00014710126801485444, 8.081430028285006e-05, 1.0130649785457274e-05]



def fitpolyexp(densities, binning):
  binning = numpy.array(binning)

  def polyexp(xs, a, b, c):
    coeffs = [a, b, c]
    return numpy.exp(numpy.polyval(coeffs, xs))

  centers = (binning[1:] + binning[:-1]) / 2

  optcoeff = optimize.curve_fit(polyexp, centers, densities, p0=numpy.zeros((3,)))[0]

  def f(xs):
    return polyexp(xs, optcoeff[0], optcoeff[1], optcoeff[2])

  return f



MCfitlow = interpolate.interp1d(centers, toppthist, kind="linear", fill_value=toppthist[0], bounds_error=False)
datafitlow = interpolate.interp1d(centers, toppt["densities"], kind="linear", fill_value=toppt["densities"][0], bounds_error=False)


MCfithigh = fitpolyexp(toppthist[BINSPLIT:], ptbinning[BINSPLIT:])
datafithigh = fitpolyexp(toppt["densities"][BINSPLIT:], ptbinning[BINSPLIT:])


def MCfit(xs):
  return numpy.where(xs < ptbinning[BINSPLIT], MCfitlow(xs), MCfithigh(xs))

def datafit(xs):
  return numpy.where(xs < ptbinning[BINSPLIT], datafitlow(xs), datafithigh(xs))



if __name__ == "__main__":
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
  fig.savefig("dataMC.pdf")

