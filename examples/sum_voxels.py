import matplotlib
import matplotlib.pylab as plt
import numpy as np
import sys

################################################################################
# Open the input file
################################################################################
infilenames = sys.argv[1:]

ngals = None
nbins = None
histvals = [None,None,None]

totpts = None

infile = []
for fcount,f in enumerate(infilenames):
    infile = open(f,'r')

    vals = np.array(infile.read().split()).astype('str')

    ngals = vals[0:3].astype('int')
    nbins = vals[3:6].astype('int')
    for i in xrange(3):
        histvals[i] = vals[6+(i*3):9+(i*3)].astype('float')

    pts = vals[15:].astype('float')

    if fcount==0:
        totpts = np.zeros_like(pts)

    totpts += pts

    #pts = pts.reshape((nbins[0],nbins[1],nbins[2]))

#print totpts
totpts = totpts.reshape((nbins[0],nbins[1],nbins[2]))
#for t in totpts:
    #print t

output = ""
for n in ngals:
    output += "%d\n" % n
output += "%d %d %d\n" % (nbins[0],nbins[1],nbins[2])
for h in histvals:
    output += "%-6.3f %-6.3f %-6.3f\n" % (h[0],h[1],h[2])
for i in xrange(nbins[0]):
    for j in xrange(nbins[1]):
        for k in xrange(nbins[2]):
            #index = (nbins[0]*nbins[1])*i + nbins[1]*j + k
            #output += "%d" % (totpts[index])
            output += "%d " % (totpts[i][j][k])
        output += "\n"

print output


