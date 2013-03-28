import matplotlib
import matplotlib.pylab as plt
import numpy as np
import sys

################################################################################
# Open the input file
################################################################################
infilename = []
#infilename.append('DDD_evenbinning_CPU_10k.dat')
#infilename.append('DDR_evenbinning_CPU_10k.dat')
#infilename.append('DRR_evenbinning_CPU_10k.dat')
#infilename.append('RRR_evenbinning_CPU_10k.dat')

infilename.append('DDD_evenbinning_CPU_128bins_10k.dat')
infilename.append('DDR_evenbinning_CPU_128bins_10k.dat')
infilename.append('DRR_evenbinning_CPU_128bins_10k.dat')
infilename.append('RRR_evenbinning_CPU_128bins_10k.dat')

ddd = None
ddr = None
drr = None
rrr = None

infile = []
for fcount,f in enumerate(infilename):
    infile = open(f,'r')

################################################################################

    nbins =  [0, 0, 0]
    pts = None

    i = 0
    j = 0
    k = 0

    vals = np.array(infile.read().split()).astype('longlong')
    nbins = [vals[0]+2,vals[1]+2,vals[2]+2]
    pts = np.zeros((nbins[0],nbins[1],nbins[2]))
    print "here"
    print nbins
    print vals

    for i in range(0,nbins[0]):
        for j in range(0,nbins[1]):
            for k in range(0,nbins[2]):
                #print vals[3 + i*1 + i*nbins[1]*nbins[2] + j*nbins[1] + k]
                index = 3 + (i+1) + i*nbins[1]*nbins[2] + j*nbins[2] + k
                #print index
                pts[i][j][k] = float(vals[3 + (i+1) + i*nbins[1]*nbins[2] + j*nbins[1] + k])

    if fcount==0:
        ddd = pts.copy()
    elif fcount==1:
        ddr = pts.copy()
    elif fcount==2:
        drr = pts.copy()
    elif fcount==3:
        rrr = pts.copy()

print "Sums: -----------"
ddd_norm =  float(sum(sum(sum(ddd))))
ddr_norm =  float(sum(sum(sum(ddr))))
drr_norm =  float(sum(sum(sum(drr))))
rrr_norm =  float(sum(sum(sum(rrr))))

index0 = ddd>1000
index1 = ddr>1000
index2 = drr>1000
index3 = rrr>1000

ddd /= ddd_norm
ddr /= ddr_norm
drr /= drr_norm
rrr /= rrr_norm

good_index = index0*index1*index2*index3

'''
for x0,x1,x2,x3 in zip(ddd[good_index==0] , ddr[good_index==0] , drr[good_index==0] , rrr[good_index==0]):
    print x0,x1,x2,x3
'''

tpcf = np.zeros((nbins[0],nbins[1],nbins[2]))
tpcf[good_index] = (ddd[good_index] - (3*ddr[good_index]) + (3*drr[good_index]) - rrr[good_index])/rrr[good_index]

tpcf[tpcf!=tpcf] = 0
tpcf[tpcf==np.inf] = 0

fig = []
axes = []
#tag = infilename.split('/')[-1].split('.dat')[0]
for i in range(nbins[0]-2):
    #axes.append(figure.add_subplot(nbins[0]-2,nbins[1]-2,i*(nbins[0]-2) + j + 1))
    fig.append(plt.figure())
    axes.append(fig[i].add_subplot(1,1,1))
    extent = [0,nbins[1],0,nbins[2]]
    #print extent
    #print tpcf[i+1]
    norm = matplotlib.colors.Normalize(vmin=-1.0,vmax=1.0)
    cs = axes[i].imshow(tpcf[i+1],extent=extent,interpolation='nearest',origin='lower',cmap=plt.cm.coolwarm,axes=axes[i],aspect='auto',norm=norm)
    #axes[i].imshow(vals[i],interpolation='nearest')

    #axes[i].set_zlim(-1,1)
    plt.colorbar(cs)

    name = "Plots/tpcf_fig%03d.png" % (i)
    fig[i].savefig(name)

print tpcf
'''
for t,t0,t1,t2,t3 in zip(tpcf,ddd,ddr,drr,rrr):
    print " ---------- "
    print t
    print t0
    print t1
    print t2
    print t3
'''

print "Sums: -----------"
print ddd_norm
print ddr_norm
print drr_norm
print rrr_norm


#plt.show()
