import numpy as np
import matplotlib.pylab as plt

import sys

infile = open(sys.argv[1],'r')

vals = np.array(infile.read().split()).astype('float')

ngals = int(vals[0])

ncols = 3
index = np.arange(1,ngals*ncols,ncols)
x = vals[index]
y = vals[index+1]
x = vals[index+2]

plt.plot(x,y,'o',markersize=1)
plt.show()
