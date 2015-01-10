import matplotlib
import matplotlib.pylab as plt
import numpy as np
import sys


################################################################################
def sxx2cf(s,xpts,ypts):

    cf = None
    if type(s)==np.ndarray:
        cf = np.zeros(len(s))

        for i,spt in enumerate(s):
            diff = np.abs(spt-xpts)
            m = min(diff)
            index = diff.tolist().index(m)

            cf[i] = ypts[index]
    else:
        diff = np.abs(s-xpts)
        m = min(diff)
        index = diff.tolist().index(m)

        cf = ypts[index]

    return cf

################################################################################
def bin2val(ibin,lo,hi,nbins):

    width = (hi-lo)/float(nbins)

    val = lo + (width*ibin) + (width/2.0)

    return val

################################################################################

################################################################################
def sqstheta2s1s2s3(s,qs,theta):

    s1 = s
    s2 = s*qs
    s3 = s*np.sqrt(1 + qs*qs - 2*qs*np.cos(theta))

    return s1,s2,s3


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
        ddd_norm =  float(ngals[0]*(ngals[1]-1)*(ngals[2]-2))/6.
    elif fcount==1:
        ddr = pts.copy()
        ddr_norm =  float(ngals[0]*(ngals[1]-1)*ngals[2])/2.
    elif fcount==2:
        drr = pts.copy()
        drr_norm =  float(ngals[0]*ngals[1]*(ngals[2]-1))/2.
    elif fcount==3:
        rrr = pts.copy()
        rrr_norm =  float(ngals[0]*(ngals[1]-1)*(ngals[2]-2))/6.

snbins = nbins[0]
slo = histvals[0][0]
shi = histvals[0][1]
swidth = histvals[0][2]

qsnbins = nbins[1]
qslo = histvals[1][0]
qshi = histvals[1][1]
qswidth = histvals[1][2]

thetanbins = nbins[2]
thetalo = histvals[2][0]
thetahi = histvals[2][1]
thetawidth = histvals[2][2]

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

tpcf = np.zeros((snbins,qsnbins,thetanbins))

#tpcf[good_index] = (ddd[good_index] - (3*ddr[good_index]) + (3*drr[good_index]) - rrr[good_index])/rrr[good_index]
tpcf = (ddd - (3*ddr) + (3*drr) - rrr)/rrr

tpcf[tpcf!=tpcf] = 0
tpcf[tpcf==np.inf] = 0

#print tpcf
#print tpcf[tpcf>0]

figall = plt.figure()
figtpcf = plt.figure()
figredtpcf = plt.figure()

ax0 = figall.add_subplot(4,1,1)
ax1 = figall.add_subplot(4,1,2)
ax2 = figall.add_subplot(4,1,3)
ax3 = figall.add_subplot(4,1,4)

title = "s=%4.1f-%4.1f Mpc" % (sbin*swidth+slo,(sbin+1)*swidth+slo)
axtpcf = figtpcf.add_subplot(1,1,1)
axtpcf.set_title(title)

title = "Reduced s=%4.1f-%4.1f Mpc" % (sbin*swidth+slo,(sbin+1)*swidth+slo)
axredtpcf = figredtpcf.add_subplot(1,1,1)
axredtpcf.set_title(title)

################################################################################
# Read in 2pcf stuff.
infile = open("master_2pcf_reference.dat")
vals = np.array(infile.read().split()).astype('float')
nvals = len(vals)
index = np.arange(0,nvals,2)
x2pcf = vals[index]
y2pcf = vals[index+1]
print x2pcf
print y2pcf
#exit()
################################################################################

x = np.arange(thetalo,thetahi,thetawidth)

print x

sval = bin2val(sbin,slo,shi,snbins)

for qs in qsbin:

    qsval = bin2val(qs,qslo,qshi,qsnbins)
    print "qsval: %f" % (qsval)

    ax0.plot(x,ddd[sbin][qs],'o')
    ax1.plot(x,ddr[sbin][qs],'o')
    ax2.plot(x,drr[sbin][qs],'o')
    ax3.plot(x,rrr[sbin][qs],'o')

    thetaval = np.pi*tpcf[sbin][qs]

    s12,s23,s31 = sqstheta2s1s2s3(sval,qsval,thetaval)
    #print s12,s23,s31
    s12cf = sxx2cf(s12,x2pcf,y2pcf)
    s23cf = sxx2cf(s23,x2pcf,y2pcf)
    s31cf = sxx2cf(s31,x2pcf,y2pcf)
    denominator = s12cf*s23cf + s23cf*s31cf + s31cf*s12cf

    label = "qs=%3.1f-%3.1f" % (qs*qswidth+qslo,(qs+1)*qswidth+qslo)
    axtpcf.plot(x,tpcf[sbin][qs],'o-',label=label)
    plt.xlabel(r'$\theta/\pi$',fontsize=24)
    plt.xlim(0,1.5)

    label = "qs=%3.1f-%3.1f" % (qs*qswidth+qslo,(qs+1)*qswidth+qslo)
    axredtpcf.plot(x,tpcf[sbin][qs]/denominator,'o-',label=label)
    print label
    print denominator
    print tpcf[sbin][qs]/denominator
    plt.xlabel(r'$\theta/\pi$',fontsize=24)
    plt.xlim(0,1.5)

plt.legend()
axtpcf.legend()
plt.tight_layout()

print "Sums: -----------"
print ddd_norm
print ddr_norm
print drr_norm
print rrr_norm

plt.show()
