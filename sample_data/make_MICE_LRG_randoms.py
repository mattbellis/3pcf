## I'm going to go ahead and assume there are no survey masks involved here.
## so, I'm just going to randomly scatter 'rans=doms' in x,y,z.
## note that we only actually use x/y/z in our code, so ra/dec/z are irrelevant. 

import copy
from ROOT import *
import time, math
import numpy as np
from array import array

nbins = 10
ra, dec, z = [], [], []
x, y, z = [], [], []


### make the log x axis.
xmin = 0.01
xmax = 250
logxmin = math.log10(xmin)
logxmax = math.log10(xmax)
binwidth = (logxmax-logxmin)/nbins
xbins = []
ybins = []
for i in range(nbins+1):
    xbins.append( xmin + math.pow(10, (logxmin+i) * binwidth))
    ybins.append( xmin + math.pow(10, (logxmin+i) * binwidth))

print nbins
print xbins
xbins = array('f',xbins)
ybins = array('f',ybins)

h_x = TH1F('','',nbins, 0, 250)
h_y = TH1F('','',nbins, 400, 1500)
h_z = TH1F('','',nbins, 0, 250)

h_xy = TH2F('','', nbins, 0, 500, nbins, 400, 1500)
h_xz = TH2F('','', nbins, 0, 500, nbins, 0, 500)
h_yz = TH2F('','',nbins,400,1500, nbins,0,500)

h_xyz = TH3F('','', nbins,0,500, nbins,400,1500, nbins,0,500)

ngals = 0
for line in open('MICE_LRGs_20degx20deg.dat'):
    cols = line.split()
    ngals+=1
    
    x.append(float(cols[4]))
    y.append(float(cols[5]))
    z.append(float(cols[6]))
    h_x.Fill(float(cols[4]))
    h_y.Fill(float(cols[5]))
    h_z.Fill(float(cols[6]))

    h_xy.Fill(float(cols[4]), float(cols[5]))
    h_xz.Fill(float(cols[4]), float(cols[6]))
    h_yz.Fill(float(cols[5]), float(cols[6]))

    h_xyz.Fill(float(cols[4]), float(cols[5]), float(cols[6]))

xx = copy.deepcopy(x)
xx.sort()
yy = copy.deepcopy(y)
yy.sort()
zz = copy.deepcopy(z)
zz.sort()
print xx[0], xx[int(len(xx)-1)]
print yy[0], yy[int(len(yy)-1)]
print zz[0], zz[int(len(zz)-1)]


h_xyz.Draw()


c1 = TCanvas()
h_xy.Draw('colz')
c2 = TCanvas()
h_xz.Draw('colz')
c3 = TCanvas()
h_yz.Draw('colz')

             
    
    
### now get the random distribution
## gonna do 16x randoms

## there'll be some binning artifacts, cos there's so few original gals to work from. That's a shame!

hh_xyz = TH3F('', '', 50,0,250, 50,400,1500, 50,0,250)

flatfile = open("flat_MICE_LRGs_20degx20deg.dat","w")
for i in range(ngals*16):
    
    rx, ry, rz = Double(0.), Double(0.), Double(0.)
    h_xyz.GetRandom3(rx, ry, rz)
    ## id, ra, dec, z are irrelevant
    print >> flatfile, "00", "0","0", "0", rx, ry, rz
    hh_xyz.Fill(rx, ry, rz)
    
flatfile.close()
print rx, ry, rz


ccc = TCanvas()
hh_xyz.Draw()

time.sleep(300)
