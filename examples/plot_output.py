import matplotlib.pylab as plt
import numpy as np
import sys

################################################################################
# Open the input file
################################################################################
infilename = None
if sys.argv[1] is not None:
    infilename = sys.argv[1]
else:
    print "Need to pass in an input file on the command-line."
infile = open(infilename,"r")

################################################################################

nbins =  [0, 0, 0]
pts = None

i = 0
j = 0
k = 0

vals = np.array(infile.read().split()).astype('int')
nbins = [vals[0]+2,vals[1]+2,vals[2]+2]
pts = np.zeros((nbins[0],nbins[1],nbins[2]))
print nbins
print vals

for i in range(0,nbins[0]):
    for j in range(0,nbins[1]):
        for k in range(0,nbins[2]):
            #print vals[3 + i*1 + i*nbins[1]*nbins[2] + j*nbins[1] + k]
            index = 3 + (i+1) + i*nbins[1]*nbins[2] + j*nbins[2] + k
            #print index
            pts[i][j][k] = vals[3 + (i+1) + i*nbins[1]*nbins[2] + j*nbins[1] + k]



print pts

fig = []
axes = []
tag = infilename.split('/')[-1].split('.dat')[0]
for i in range(nbins[0]-2):
    #axes.append(figure.add_subplot(nbins[0]-2,nbins[1]-2,i*(nbins[0]-2) + j + 1))
    fig.append(plt.figure())
    axes.append(fig[i].add_subplot(1,1,1))
    extent = [0,nbins[1],0,nbins[2]]
    print extent
    print pts[i+1]
    axes[i].imshow(pts[i+1],extent=extent,interpolation='nearest',origin='lower',cmap=plt.cm.coolwarm,axes=axes[i],aspect='auto')
    #axes[i].imshow(vals[i],interpolation='nearest')

    name = "Plots/%sfig%03d.png" % (tag,i)
    fig[i].savefig(name)


#plt.show()
