import numpy as np
import matplotlib.pylab as plt
from mpl_toolkits.mplot3d import Axes3D
import scipy.spatial
from time import time

################################################################################
# Use this to define the ranges and boundaries for the volumetric 
# ``chunks" of galaxies.
################################################################################
def define_ranges(vec_ranges,maxsep=150):

    ndim = len(vec_ranges)
    #print ndim

    vec_nbins = []
    vec_binwidths = []
    for i in range(ndim):

        r = vec_ranges[i][1]-vec_ranges[i][0]
        nbins = int(r/maxsep)# - 1
        binwidth = r/nbins

        vec_nbins.append(nbins)
        vec_binwidths.append(binwidth)

        #print "------------"
        #print nbins
        #print "----"
        '''
        for j in range(0,nbins+1):
            print binwidth*j
        '''

    return vec_nbins,vec_binwidths

################################################################################
# Start the tests
################################################################################
ngals = 100000
maxsep = 300
total_calcs = float(ngals*ngals - ngals)/2.

print "Running with %d galaxies\n" % (ngals)

# Range of SDSS DR10
x = (-2000.,140.)
y = (-1720.,1780.)
z = (755.,2500.)

################################################################################
# Generate some fake galaxies
################################################################################
print "Generating the galaxies...."
start = time()
np.random.seed(1)

gal0 = np.random.random((3,ngals))
gal0[0] = (x[1]-x[0])*gal0[0] + x[0]
gal0[1] = (y[1]-y[0])*gal0[1] + y[0]
gal0[2] = (z[1]-z[0])*gal0[2] + z[0]

gal0 = gal0.transpose()
print "Generated and transposed the galaxies...."
print "time: %f" % (time()-start)
#print gal0

print "\n"

################################################################################
# Do the calculations (DD for example) the normal way.
################################################################################
print "-------------------------------------\n"
'''
print "Doing the full calculations.\n"
print "Calculating the distances."
start = time()
starta = time()
# Do it the normal way.
#dists = scipy.spatial.distance.cdist(gal0,gal0)
dists = scipy.spatial.distance.pdist(gal0)
print "Calculated the distances.\ttime: %f" % (time()-starta)

print "Cleaning up the distances...."
starta = time()
dists = dists.flatten()
index = dists<maxsep
index *= dists!=0
cleaned_up_distances = np.unique(dists[index])
print "Cleaned up distances.\t\ttime: %f\n" % (time()-starta)

total_calcs = len(dists)
histogrammed_calcs = len(cleaned_up_distances)
print "Total \t\t\t\ttime: %f\n" % (time()-start)
print "Total/histogrammed calcs: %d %d %f\n" % (total_calcs,histogrammed_calcs,histogrammed_calcs/float(total_calcs))
'''

#exit()


################################################################################
# Do it broken up into smaller chunks
################################################################################
print "-------------------------------------\n"
print "Doing it the chunked way.\n"
print "Indexing the galaxies......"
start = time()
starta = time()

# Get the bins and bin widths
vec_ranges = [x,y,z]
vec_nbins,vec_binwidths = define_ranges(vec_ranges,maxsep=maxsep)
# Total number of voxels.
totbins = 1
for v in vec_nbins:
    totbins *= v

print "--------"
print "totbins:"
print vec_nbins,totbins
print "--------"

# Calculate the voxel location (3 numbers) for each galaxy.
allbins = []
for gal in gal0:
    bins = []
    for i,g in enumerate(gal):
        bin = int((g-vec_ranges[i][0])/vec_binwidths[i])
        bins.append(bin)
    #print gal,bins,vec_ranges,vec_binwidths #vec_nbins
    allbins.append(np.array(bins))

allbins = np.array(allbins)
#print allbins
allbins = allbins.transpose()
#print allbins
print "Indexed the galaxies.\t\ttime: %f\n" % (time()-starta)

#print gal0
#print gal0[allbins[0]==0]
#print len(gal0[allbins[0]==0])

print "Divide up the galaxies and copy over into arrays for easier manipulations....."
chunked_dists = []
#chunked_dists = np.array([])
totgals = 0
chunked_galaxies = []

# First divide up the galaxies
print "Dividing up the galaxies...."
starta = time()
for i in range(vec_nbins[0]):
    chunked_galaxies.append([])
    for j in range(vec_nbins[1]):
        chunked_galaxies[i].append([])
        for k in range(vec_nbins[2]):
            #chunked_galaxies[i][j].append([])

            #print i,j,k

            index = allbins[0]==i
            index *= allbins[1]==j
            index *= allbins[2]==k

            # Primary galaxies
            galaxies0 = gal0[index]
            #print "here"
            #print type(galaxies0)
            chunked_galaxies[i][j].append(galaxies0)

            totgals += len(galaxies0)
print "Divided the galaxies.\t\ttime: %f\n" % (time()-starta)

# Then do the calculations
print "Calculating the distances......."
starta = time()
for i in range(vec_nbins[0]):
    for j in range(vec_nbins[1]):
        for k in range(vec_nbins[2]):

            galaxies0 = chunked_galaxies[i][j][k]

            #print galaxies0
            #print type(galaxies0)

            # Secondary galaxies
            galaxies1 = []
            #'''
            index = np.zeros(ngals).astype('bool')
            for ii in range(-1,2):
                for jj in range(-1,2):
                    for kk in range(-1,2):
                        #index1 = allbins[0]==i+ii
                        #index1 *= allbins[1]==j+jj
                        #index1 *= allbins[2]==k+kk

                        #print i,i+ii,j,j+jj,k,k+kk
                        #index += index1
                        iindex = i+ii
                        jindex = j+jj
                        kindex = k+kk
                        if iindex<vec_nbins[0] and jindex<vec_nbins[1] and kindex<vec_nbins[2]:
                            galaxies1 += chunked_galaxies[iindex][jindex][kindex].tolist()

            #galaxies1 = gal0[index]
            #print galaxies1
            galaxies1 = np.array(galaxies1)
            #print galaxies1

            #print len(galaxies0),len(galaxies1)
            #print i,j,k,galaxies0,galaxies1
            if len(galaxies0)>0 and len(galaxies1)>0:
                dists = scipy.spatial.distance.cdist(galaxies0,galaxies1)
                #print dists
                dists = dists.flatten()
                #print len(dists)
                #print "----"
                #print len(dists)
                index = dists<maxsep
                index*= dists!=0
                dists = dists[index]
                #print len(dists)
                chunked_dists += dists.tolist()
                #print len(chunked_dists)
                #print len(chunked_dists)
                #chunked_dists = np.append(chunked_dists,dists)
            #'''

    #print len(chunked_dists)

            #print totgals

print "Calculated the distances.\ttime: %f\n" % (time()-starta)

print "Cleaning up the distances...."
starta = time()
chunked_dists = np.array(chunked_dists)
#chunked_dists = chunked_dists.flatten()
index = chunked_dists<maxsep
index *= chunked_dists!=0
cleaned_up_distances = np.unique(chunked_dists[index])
print "Cleaned up distances.\t\ttime: %f\n" % (time()-starta)

total_calcs = len(chunked_dists)
histogrammed_calcs = len(cleaned_up_distances)
print "Total \t\t\t\ttime: %f\n" % (time()-start)
print "Total/histogrammed calcs: %d %d %f\n" % (total_calcs,histogrammed_calcs,histogrammed_calcs/float(total_calcs))

#exit()

#vec_nbins,vec_binwidths = define_ranges(vec_ranges)


#gal0 = gal0.transpose()
#fig = plt.figure()
#ax = fig.add_subplot(111, projection='3d')
#ax.scatter(gal0[0], gal0[1], gal0[2], s=1.0)

#plt.hist(chunked_dists,bins=100)
plt.show()
