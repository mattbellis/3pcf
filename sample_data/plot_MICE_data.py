import numpy as np
import matplotlib.pylab as plt

import sys

infile = open(sys.argv[1],'r')

vals = np.array(infile.read().split()).astype('float')

ncols = 7
ngals = len(vals)/ncols

print ngals

index = np.arange(0,(ngals*ncols)-1,ncols)
#print index

id = vals[index].astype('int')
ra = vals[index+1]
dec = vals[index+2]
zredshift = vals[index+3]
x = vals[index+4]
y = vals[index+5]
z = vals[index+6]

print "ra        - min/max: %f %f" % (min(ra),max(ra))
print "dec       - min/max: %f %f" % (min(dec),max(dec))
print "zredshift - min/max: %f %f" % (min(zredshift),max(zredshift))
print "x         - min/max: %f %f" % (min(x),max(x))
print "y         - min/max: %f %f" % (min(y),max(y))
print "z         - min/max: %f %f" % (min(z),max(z))
r = np.sqrt(x*x + y*y + z*z)
print "r         - min/max: %f %f" % (min(r),max(r))

################################################################################

markersize=0.5
alpha=0.5

plt.figure(figsize=(12,4))

plt.subplot(1,3,1)
plt.plot(x,y,'o',markersize=markersize,alpha=alpha)
plt.xlabel('x')
plt.ylabel('y')

plt.subplot(1,3,2)
plt.plot(x,z,'o',markersize=markersize,alpha=alpha)
plt.xlabel('x')
plt.ylabel('z')

plt.subplot(1,3,3)
plt.plot(y,z,'o',markersize=markersize,alpha=alpha)
plt.xlabel('y')
plt.ylabel('z')

plt.tight_layout()


################################################################################

plt.figure(figsize=(12,8))

plt.subplot(2,3,1)
plt.plot(ra,dec,'o',markersize=markersize,alpha=alpha)
plt.xlabel('ra')
plt.ylabel('dec')

plt.subplot(2,3,2)
plt.plot(ra,zredshift,'o',markersize=markersize,alpha=alpha)
plt.xlabel('ra')
plt.ylabel('zredshift')

plt.subplot(2,3,3)
plt.plot(dec,zredshift,'o',markersize=markersize,alpha=alpha)
plt.xlabel('dec')
plt.ylabel('zredshift')

################################

plt.subplot(2,3,4)
plt.plot(ra,np.cos(np.deg2rad(90-dec)),'o',markersize=markersize,alpha=alpha)
plt.xlabel('ra')
plt.ylabel(r'$\cos(\pi/2-dec)$')

plt.subplot(2,3,5)
plt.plot(ra,zredshift,'o',markersize=markersize,alpha=alpha)
plt.xlabel('ra')
plt.ylabel('zredshift')

plt.subplot(2,3,6)
plt.plot(np.cos(np.deg2rad(90-dec)),zredshift,'o',markersize=markersize,alpha=alpha)
plt.xlabel(r'$\cos(\pi/2-dec)$')
plt.ylabel('zredshift')


plt.tight_layout()


plt.show()
