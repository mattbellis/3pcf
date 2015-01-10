import matplotlib
import matplotlib.pylab as plt
import numpy as np
import sys

################################################################################
# Open the input file
################################################################################
infilename = [None,None,None]

infilename[0] = sys.argv[1]
infilename[1] = infilename[0].replace('DD','DR')
infilename[2] = infilename[0].replace('DD','RR')

dd = None
dr = None
rr = None

dd_norm = None
dr_norm = None
rr_norm = None

ngals = None
nbins = None
histvals = [None,None,None]

infile = []
for fcount,f in enumerate(infilename):
    print f
    infile = open(f,'r')

################################################################################

    vals = np.array(infile.read().split()).astype('str')

    ngals = vals[0:2].astype('int')
    nbins = vals[2].astype('int')
    histvals = vals[3:6].astype('float')

    pts = vals[6:].astype('float')

    if fcount==0:
        dd = pts.copy()
        dd_norm =  float(ngals[0]*ngals[0]-ngals[0])/2.
    elif fcount==1:
        dr = pts.copy()
        dr_norm =  float(ngals[0]*ngals[1])
    elif fcount==2:
        rr = pts.copy()
        rr_norm =  float(ngals[0]*ngals[0]-ngals[0])/2.


# Only wory about entries greater than 1000
'''
index0 = ddd>1000
index1 = ddr>1000
index2 = drr>1000
index3 = rrr>1000

good_index = index0*index1*index2*index3

'''

dd /= dd_norm
dr /= dr_norm
rr /= rr_norm

#print ddd
#print ddd_norm

#tpcf[good_index] = (ddd[good_index] - (3*ddr[good_index]) + (3*drr[good_index]) - rrr[good_index])/rrr[good_index]
tpcf = (dd - (2*dr) + rr)/rr

tpcf[tpcf!=tpcf] = 0
tpcf[tpcf==np.inf] = 0

#print tpcf
#print tpcf[tpcf>0]

figall = plt.figure()
figtpcf = plt.figure()

ax0 = figall.add_subplot(3,1,1)
ax1 = figall.add_subplot(3,1,2)
ax2 = figall.add_subplot(3,1,3)

#lo = histvals[0][0]
#hi = histvals[0][1]
#width = histvals[0][2]
title = "2pcf"

axtpcf = figtpcf.add_subplot(1,1,1)
axtpcf.set_title(title)

lo = histvals[0]
hi = histvals[1]
width = histvals[2]
print lo,hi,width
x = np.arange(lo,hi,width)

print x
print len(x)
print len(dd)

ax0.plot(x,dd,'o')
ax1.plot(x,dr,'o')
ax2.plot(x,rr,'o')


label = "plot"
axtpcf.plot(x,tpcf,'o-',label=label)

plt.xlabel(r's (Mpc)',fontsize=24)
#plt.xlim(0,1.5)
plt.yscale('log')

plt.legend()
plt.tight_layout()

print "Sums: -----------"
print dd_norm
print dr_norm
print rr_norm

plt.show()
