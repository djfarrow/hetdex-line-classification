"""

Generate a data cube containing the equivalent 
width distributions as a function of wavelength. For
use in the new "empirical" model of EW.

AUTHOR: Daniel Farrow (MPE)

"""
from __future__ import print_function
import argparse
import matplotlib.pyplot as plt
from ConfigParser import RawConfigParser
from numpy import( linspace, histogram, exp, digitize, mean, polyfit, polyval, zeros, log10, power, 
                   trapz )
from scipy.optimize import leastsq
from astropy.table import Table
from astropy.io import fits
from line_classifier.lfs_ews.equivalent_width import EquivalentWidthAssigner

def diff_func(w0_vec, ew_bins, ew_vals):
    """ Tool for fitting models to the distribution """
    w0 = w0_vec[0]
    model = (1.0/w0)*exp(-1.0*ew_bins/w0)
    
    return ew_vals - model 

   
parser = argparse.ArgumentParser(description="Fit observed EW distribution from mock. Return new w0 values")
parser.add_argument("config", help="Config filename (same config as rest of code)")
parser.add_argument("cats", nargs="*",  help="Catalogue filename")
parser.add_argument("out", help="Output the EW distributions here")
parser.add_argument("--oii", action='store_true', help="An OII catalogue?")
parser.add_argument("--plot", action='store_true', help="Plot the distributions")
opts = parser.parse_args()

config = RawConfigParser()
config.read(opts.config)

ew_max = config.getfloat("InterpolatedEW", "oii_ew_max") 

if opts.oii :
    newbins = 100
    ew = EquivalentWidthAssigner.from_config(config, 'OII_EW')
    ew_bins = linspace(log10(0.5), log10(ew_max), newbins)
else:
    newbins = 100
    ew = EquivalentWidthAssigner.from_config(config, 'LAE_EW')
    ew_bins = linspace(log10(0.5), log10(ew_max), newbins)

nbins = 20
# Redshift range over which OII emitters 
# are actually detectable
if opts.oii:
    # Bins smaller for OII. ### Want this?
    bins = linspace(0.05, 0.5, nbins)
else:
    bins = linspace(2.18, 3.51, nbins)

zbcens = 0.5*(bins[:-1] + bins[1:])

lew_bcens = 0.5*(ew_bins[:-1] + ew_bins[1:])
ew_bcens = power(10, lew_bcens)
ew_bsize = power(10, ew_bins[1:]) - power(10, ew_bins[0:-1])

hist_cube = zeros((nbins - 1, newbins - 1))

# Remember number above EW range
nabove = zeros(nbins - 1)

# Only need this for the plot
if opts.plot:
    hist_cube_true = zeros((nbins - 1, newbins - 1))

for fname in opts.cats:

    print(fname)

    table = Table.read(fname)
    
    ew_obs = table["EW_OBS"]
    ew_true = table["EW_TRUE"]
    z = table["z_redshift_space_true"]
    ew_obs_rest = ew_obs/(1.0 + z)
    ew_true_rest = ew_true/(1.0 + z)
    
    indices = digitize(z, bins=bins)
    
    # In redshift bins fit an exponential
    for index in range(1, nbins):
    
        tew_obs_rest = ew_obs_rest[(indices == index) & (ew_obs_rest > 0)]
        hist_obs = histogram(log10(tew_obs_rest), bins=ew_bins)[0]
    
        hist_cube[index - 1, :] += hist_obs
        nabove[index - 1] += len(tew_obs_rest[tew_obs_rest > ew_max])

        # Store this for a plot later
        if opts.plot:    
            tew_true_rest = ew_true_rest[indices == index]
            hist_true = histogram(log10(tew_true_rest), bins=ew_bins)[0]
            hist_cube_true[index - 1, :] += hist_true


for index, na in enumerate(nabove, 1):

    # Remember to include sources outside of range in normalisation
    norm_obs = (sum(hist_cube[index - 1, :]) + na)*ew_bsize

    print(sum(hist_cube[index - 1, :]), na)

    hist_cube[index - 1, :] /= norm_obs

    print(trapz(hist_cube[index - 1, :], x=ew_bcens))

    # Avoid NaNs like the plague
    hist_cube[index - 1, norm_obs < 1e-99] = 0
 


# Plot the EW distributions
if opts.plot:   

    for index in range(2, nbins, 2):
 
        norm_true = sum(hist_cube_true[index - 1, :])*ew_bsize
        hist_cube_true[index - 1,:] /= norm_true
  
        # Base these of the last catalogue used (is ther a better approach??)
        meanz = mean(z[indices == index])
        ew_func = ew.ew_function(meanz)  

        # Plot of log10 scale
        trans = lambda x : log10(x)

        ax = plt.subplot(3, 3, index/2)
        ax.plot(ew_bcens, trans(hist_cube[index - 1, :]), 'k-', marker="*", label="Observed EW") 
        ax.plot(ew_bcens, trans(hist_cube_true[index - 1, :]), 'r-', marker="*", label="True EW")
        ax.plot(ew_bcens, trans(ew_func(ew_bcens)), 'r:', linewidth=3.0, label="Input EW")

        w0fit_true, info = leastsq(diff_func, [50.0], args=(ew_bcens, hist_cube_true[index - 1, :]))
        w0fit, info = leastsq(diff_func, [50.0], args=(ew_bcens, hist_cube[index - 1, :]))
        ax.plot(ew_bcens, trans((1.0/w0fit[0])*exp(-1.0*ew_bcens/w0fit[0])), 'k:', linewidth=3.0)     

        plt.text(10.0, ax.get_ylim()[0]*0.9, "{:2.1f}<z<{:2.1f}".format(bins[index - 1], bins[index]), fontsize=14.0)

        if index/2 == 1 or index/2 == 4 or index/2 == 7:
            ax.set_ylabel('log10(dN/dEW)')
        ax.set_xlabel('Equivalent Width')


        if index/2 == 3:
            plt.legend(frameon=False, loc="upper right")

    plt.show()

print("Number of sources above max EW in z-bins: ", nabove)

hdu1 = fits.PrimaryHDU(hist_cube)
hdu2 = fits.ImageHDU(data=ew_bcens, name="EW_BCENS")
hdu3 = fits.ImageHDU(data=zbcens, name="REDSHIFT")

hdul = fits.HDUList([hdu1, hdu2, hdu3])
hdul.writeto(opts.out)


