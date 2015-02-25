import sys, random
import numpy as np

infile = sys.argv[1]

outfile = open(sys.argv[2], 'w')





## how many gals in the original file?
totgals = 0
decs, ras = [], []
for line in open(infile):
    if 'id' in line:
        continue
    totgals+=1
    cols = line.split(',')
    ras.append(float(cols[1]))
    decs.append(float(cols[2]))

    
rands = random.sample(xrange(totgals), 1000)


i = -1
for line in open(infile):
    cols = line.split(',')
    if 'id' in line:
        continue
    i+=1
    ## ra and dec are in the range 0 -> 90 degrees. 
    if float(cols[1])>20 or float(cols[2])>20:
        continue
    print >> outfile, cols[0], cols[1], cols[2], cols[3], cols[4], cols[5],  cols[6][:-1]

outfile.close()
