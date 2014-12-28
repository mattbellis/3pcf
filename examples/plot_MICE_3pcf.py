import matplotlib
import matplotlib.pylab as plt
import numpy as np
import sys

################################################################################
# Open the input file
################################################################################
infilename = [None,None,None,None]

infilename[0] = sys.argv[1]
infilename[1] = infilename[0].replace('DDD','DDR')
infilename[2] = infilename[0].replace('DDD','DRR')
infilename[3] = infilename[0].replace('DDD','RRR')

sbin = int(sys.argv[2])
qsbin = np.array(sys.argv[3:]).astype('int')

ddd = None
ddr = None
drr = None
rrr = None

ddd_norm = None
ddr_norm = None
drr_norm = None
rrr_norm = None

ngals = None
nbins = None
histvals = [None,None,None]

infile = []
for fcount,f in enumerate(infilename):
    infile = open(f,'r')

################################################################################

    vals = np.array(infile.read().split()).astype('str')

    ngals = vals[0:3].astype('int')
    nbins = vals[3:6].astype('int')
    for i in xrange(3):
        histvals[i] = vals[6+(i*3):9+(i*3)].astype('float')

    pts = vals[15:].astype('float')

    pts = pts.reshape((nbins[0],nbins[1],nbins[2]))

    if fcount==0:
        ddd = pts.copy()
        ddd_norm =  float(ngals[0]*(ngals[1]-1)*(ngals[2]-2))/6
    elif fcount==1:
        ddr = pts.copy()
        ddr_norm =  float(ngals[0]*(ngals[1]-1)*ngals[2])/2
    elif fcount==2:
        drr = pts.copy()
        drr_norm =  float(ngals[0]*ngals[1]*(ngals[2]-1))/2
    elif fcount==3:
        rrr = pts.copy()
        rrr_norm =  float(ngals[0]*(ngals[1]-1)*(ngals[2]-2))/6


# Only wory about entries greater than 1000
'''
index0 = ddd>1000
index1 = ddr>1000
index2 = drr>1000
index3 = rrr>1000

good_index = index0*index1*index2*index3

'''

print rrr

ddd /= ddd_norm
ddr /= ddr_norm
drr /= drr_norm
rrr /= rrr_norm

#print ddd
#print ddd_norm

tpcf = np.zeros((nbins[0],nbins[1],nbins[2]))

#tpcf[good_index] = (ddd[good_index] - (3*ddr[good_index]) + (3*drr[good_index]) - rrr[good_index])/rrr[good_index]
tpcf = (ddd - (3*ddr) + (3*drr) - rrr)/rrr

tpcf[tpcf!=tpcf] = 0
tpcf[tpcf==np.inf] = 0

#print tpcf
#print tpcf[tpcf>0]

figall = plt.figure()
figtpcf = plt.figure()

ax0 = figall.add_subplot(4,1,1)
ax1 = figall.add_subplot(4,1,2)
ax2 = figall.add_subplot(4,1,3)
ax3 = figall.add_subplot(4,1,4)

lo = histvals[0][0]
hi = histvals[0][1]
width = histvals[0][2]
title = "s=%4.1f-%4.1f Mpc" % (sbin*width+lo,(sbin+1)*width+lo)

axtpcf = figtpcf.add_subplot(1,1,1)
axtpcf.set_title(title)

lo = histvals[2][0]
hi = histvals[2][1]
width = histvals[2][2]
x = np.arange(lo,hi,width)

print x

for qs in qsbin:
    ax0.plot(x,ddd[sbin][qs],'o')
    ax1.plot(x,ddr[sbin][qs],'o')
    ax2.plot(x,drr[sbin][qs],'o')
    ax3.plot(x,rrr[sbin][qs],'o')


    lo = histvals[1][0]
    hi = histvals[1][1]
    width = histvals[1][2]
    label = "qs=%3.1f-%3.1f" % (qs*width+lo,(qs+1)*width+lo)
    axtpcf.plot(x,tpcf[sbin][qs],'o-',label=label)

    plt.xlabel(r'$\theta/\pi$',fontsize=24)
    plt.xlim(0,1.5)

plt.legend()
plt.tight_layout()

print "Sums: -----------"
print ddd_norm
print ddr_norm
print drr_norm
print rrr_norm

print histvals[0]
print histvals[1]
print histvals[2]

plt.show()
