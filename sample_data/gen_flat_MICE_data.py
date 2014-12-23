import os, math, sys
import numpy as np
import csv
import sys

### I want to take ra, dec, z and convert into x,y,z in cartesian coords. This is a simplified version of the conversion code in our ccogs github repository:  https://github.com/djbard/ccogs/blob/master/angular_correlation/utility_scripts/convertRADECZ_XYZ.py
# From Debbie Bard 

################################################################################
def write_output_file(id,ra,dec,zredshift,x,y,z,file_name):

    #zip(energy,time_stamps,rise_times)
    with open(file_name,'w') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerows(zip(id,ra,dec,zredshift,x,y,z))
        f.close()

################################################################################


################################################################################
def getRedshiftMpcFit():
    ### first, I'm gonna plot my cosmo-dependent redshift/distance relation
    redshift_model = [.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5]
    Mpc_model = [209.0, 413.4, 613.5, 808.8, 999.4, 1185.3, 1366.5, 1542.9, 1714.6, 1881.7]

    #gr_model = TGraph(len(redshift_model), np.array(redshift_model), np.array(Mpc_model))
    #gr_model.Fit("pol2")
    #fitresults = gr_model.GetFunction('pol2')

    #p0 = fitresults.GetParameter(0)
    #p1 = fitresults.GetParameter(1)
    #p2 = fitresults.GetParameter(2)
    #print "polynomial fit params: p0 =", p0, ", p1 =", p1, ", p2 =", p2

    params,cov = np.polyfit(redshift_model,Mpc_model,2,cov=True)
    print params
    print cov
    p0 = params[2]
    p1 = params[1]
    p2 = params[0]

    print "polynomial fit params: p0 =", p0, ", p1 =", p1, ", p2 =", p2

    return p0, p1, p2





###############################################
### convert redshift z into Mpc
def getMpcFromRedshift(redshift, p0, p1, p2):
    Mpc = p0 + p1*redshift + p2*redshift*redshift
    return Mpc


###############################################
### Convert ra/dec/Mpc coords into x/y/z coords.
def convertRaDecMpcToXYZ(ra, dec, Mpc):
    x, y, z, = 0, 0, 0
    rad = math.pi/180.0

    # Debbie
    #x = Mpc*np.sin(rad*(-1.0*dec+90))*np.cos(rad*(ra))
    #y = Mpc*np.sin(rad*(-1.0*dec+90))*np.sin(rad*(ra))
    #z = Mpc*np.cos(rad*(-1.0*dec+90))

    # MICE?
    x = Mpc*np.sin(rad*(90-dec))*np.sin(rad*(ra))
    y = Mpc*np.sin(rad*(90-dec))*np.cos(rad*(ra))
    z = Mpc*np.cos(rad*(90-dec))

    return x, y, z


#######################################
# Conversion code
#######################################
def convert(ra, dec, z):

    ### Fit 2nd order polynomial to Redshift-Mpc relation (based on standard cosmology) and get fit function params.
    # This is what Debbie used.
    #p0, p1, p2 = getRedshiftMpcFit()

    # This is what I get from the MICE stuff.
    p0,p1,p2 = -1.07199554e+00,3.01313296e+03,-6.38827998e+02

    ### convert redshift to Mpc
    zMpc = getMpcFromRedshift(z, p0, p1, p2)

    #print "zMpc: ",zMpc

    ### convert ra,dec, Mpc to x,y,z in Mpc
    xcart, ycart, zcart = convertRaDecMpcToXYZ(ra, dec, zMpc)

    #print 'converted ra=', ra, 'dec=', dec, 'and z=', z, '\n to x=', xcart, 'y=', ycart, 'z=', zcart, 'in Mpc. ' 

    return xcart, ycart, zcart


################################################################################
def main():

    ngals = int(sys.argv[1])

    id = np.arange(0,ngals,1)

    #dec = np.random.random(ngals)
    #dec = np.arccos(dec)
    #dec = 90 - np.rad2deg(dec)
    dec = (1-0.9396926)*np.random.random(ngals) # for 20x20
    dec = np.arccos(1-dec)
    dec = 20 - np.rad2deg(dec)

    #ra = 90.*np.random.random(ngals)
    ra = 20.*np.random.random(ngals) # for 20x20

    ############################################################################
    # Gen flat in the cube of z to account for volume effect.
    minz = 0.15
    maxz = 0.55
    
    minz3 = 0.15**3
    maxz3 = 0.55**3
    
    zredshift3 = (maxz3-minz3)*np.random.random(ngals) + minz3
    zredshift = zredshift3**(1./3)

    ###############################

    print "Results:"
    x,y,z = convert(ra,dec,zredshift)

    #print x,y,z
    # Give it a diffent name than ``production" files so we don't accidentlly 
    # commit different versions of a 100k line text file to git.  :)
    name = "flat_MICE_20degx20deg_%dk.dat" % (ngals/1000)
    #name = "test_flat_MICE_%dk.dat" % (ngals/1000)
    write_output_file(id,ra,dec,zredshift,x,y,z,name)


################################################################################
if __name__=="__main__":
    main()
