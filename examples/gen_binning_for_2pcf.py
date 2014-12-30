import numpy as np
s_nbins = 50
s_lo = 2.0
s_hi = 12.0

qs_nbins = 16
qs_lo = 0.9
qs_hi = 4.1

theta_nbins = 25
theta_lo = 0.
theta_hi = 1.

s_width = (s_hi-s_lo)/s_nbins
qs_width = (qs_hi-qs_lo)/qs_nbins
theta_width = (theta_hi-theta_lo)/theta_nbins

s23 = []
for s in np.arange(s_lo,s_hi,s_width):
    for qs in np.arange(qs_lo,qs_hi,qs_width):
        for theta in np.arange(theta_lo,theta_hi,theta_width):

            s12 = s
            s13 = qs*s
            s23.append(np.sqrt(s12**2 + s13**2 - 2*s12*s13*np.cos(theta)))


s23 = np.array(s23)

