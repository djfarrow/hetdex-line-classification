"""

Integrate the EW cubes to see if the 
functions normalize to 1.0 (they might not
as some sources will be outside of the EW range
considered)

AUTHOR: Daniel Farrow

"""
from __future__ import print_function

import sys
import matplotlib.pyplot as plt
from numpy import trapz
from astropy.io import fits


for fn in sys.argv[1:]:

    print(fn)
    hdus = fits.open(fn)
    integrals = [trapz(hdus[0].data[x, :], x=hdus[1].data) for x in range(hdus[2].shape[0])]
    print(integrals)
    plt.plot(hdus[1].data, hdus[0].data[7, :])


plt.show()
