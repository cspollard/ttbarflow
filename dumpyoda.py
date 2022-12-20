import yoda

yodadict = yoda.readYODA("ATLAS_2017_I1614149.yoda")

toppt = yodadict["/REF/ATLAS_2017_I1614149/d16-x01-y01"]
ttpt = yodadict["/REF/ATLAS_2017_I1614149/d20-x01-y01"]
ttm = yodadict["/REF/ATLAS_2017_I1614149/d24-x01-y01"]

def printS(s):
  binedges = []
  diffxsecs = []
  xsecs = []
  for p in s.points():
    diffxsecs.append(p.y())
    binedges.append(p.xMin())
    xsecs.append(p.y()*(p.xMax() - p.xMin()))

  binedges.append(s.points()[-1].xMax())

  print(binedges)
  print(diffxsecs)
  print(xsecs)
  return

print("toppt")
printS(toppt)
print("")

print("ttpt")
printS(ttpt)
print("")

print("ttm")
printS(ttm)
print("")
