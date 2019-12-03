"""

Generate a data cube containing the equivalent 
width distributions as a function of wavelength. For
use in the new "empirical" model of EW.

AUTHOR: Daniel Farrow (MPE)

"""
from __future__ import print_function
import argparse
from six import iteritems
from collections import OrderedDict
import matplotlib.pyplot as plt
from configparser import RawConfigParser
from numpy import(linspace, logspace, histogram, exp, digitize, mean,  
                  zeros, log10, power, meshgrid, trapz, concatenate, ndenumerate)
from astropy.table import Table
from astropy.io import fits
from line_classifier.misc.tools import generate_cosmology_from_config

   
parser = argparse.ArgumentParser(description="Generate a cube of observed flux to interpolate over")
parser.add_argument("config", help="Config filename (same config as rest of code)")
parser.add_argument("cats", nargs="*",  help="Catalogue filename")
parser.add_argument("out", help="Output the flux distributions here")
parser.add_argument("--oii", action='store_true', help="An OII catalogue?")
parser.add_argument("--plot", action='store_true', help="Plot the distributions")
opts = parser.parse_args()

config = RawConfigParser()
config.read(opts.config)
cosmo = generate_cosmology_from_config(config)

# Set up equivalent width bins
# Doesn't work well with the fast interpolation methods
fmin = 5e-18
fmax = 2e-16
f_bins = logspace(log10(fmin), log10(fmax), 50)

f_bcens = 0.5*(f_bins[:-1] + f_bins[1:])
f_bsize = f_bins[1:] - f_bins[0:-1]

nbins = 60

# Dictionary to store the bins for the histogram 
# with keys that are the parameter column name
bin_dict = OrderedDict()

# Redshift range over which OII emitters 
# are actually detectable
# the redshift bin always has to be first!
if opts.oii:
    # Binsz smaller for OII. ### Want this?
    bin_dict["z_redshift_space_true"] = linspace(0.05, 0.5, nbins)
else:
    bin_dict["z_redshift_space_true"] = linspace(2.18, 3.51, nbins)

bin_vols = cosmo.comoving_volume(bin_dict["z_redshift_space_true"][1:]).value - cosmo.comoving_volume(bin_dict["z_redshift_space_true"][:-1]).value


# Contruct a tuple with the correct number of
# indices given the number of parameters
# to bin by
shape = []
index_arrays = []
for col, bins in iteritems(bin_dict):
    shape.append(len(bins) - 1)
    index_arrays.append(range(len(bins)))
shape.append(len(f_bins) - 1)

# Have to turn it into a tuple so 
# numpy uses it as indices
shape = tuple(shape)

# Initialize an n-dimension cube to store the histograms
hist_cube = zeros(shape)

# Remember number above EW range
nabove = zeros(shape[:-1])

# Only need this for the plot
if opts.plot:
    hist_cube_true = zeros(shape)

for fname in opts.cats:

    print(fname)

    table = Table.read(fname)

    f_obs = table["FLUX_OBS"]
    f_true = table["flux"]
    z = table["z_redshift_space_true"]


    # Loop over all of the columns,
    # computing the indices in each
    # set of bins
    indices_cat = []
    for col, bins in iteritems(bin_dict):
        param = table[col]
        indices_cat.append(digitize(param, bins=bins))


    # Loop over all possible indices of output cube
    for indices, junk in ndenumerate(zeros(shape[:-1])):    

        # select sources in this n-dimensional bin
        indices_in = len(f_obs)*[True]
        for n, ind in enumerate(indices):
            indices_in = indices_in&(indices_cat[n] == (ind + 1))
            
        tf_obs = f_obs[indices_in]
        hist_obs = histogram(tf_obs, bins=f_bins)[0]

        # the redshift bins should always be first 
        hist_cube[indices][:] += hist_obs/bin_vols[indices[0]]

        # also include stuff below the range in this number
        nabove[indices] += len(tf_obs[(tf_obs > fmax)|(tf_obs < fmin)])

        # Store this for a plot later
        if opts.plot:    
            tf_true = f_true[indices_in]
            hist_true = histogram(tf_true, bins=f_bins)[0]
            hist_cube_true[indices][:] += hist_true/bin_vols[indices[0]]



# Plot the flux distributions
if opts.plot:   

    index = 1
    for indices, junk in ndenumerate(zeros(shape[:-1])):   
    
        
        labels = []
        for i, (col, bins) in enumerate(iteritems(bin_dict)):
            bin_min = bins[indices[i]]
            bin_max = bins[indices[i] + 1]
            
            if col == "z_redshift_space_true":
                zmid = 0.5*(bin_min + bin_max)

            label = "{:3.2e}<{:s}<{:3.2e}".format(bin_min, col[:1], bin_max)
            labels.append(label)
      
        # Plot of log10 scale
        trans = lambda x : log10(x)
    
        ax = plt.subplot(5, 5, index)
        #ax = plt.gca() 
        ax.plot(f_bcens, trans(hist_cube[indices][:]), 'k-', marker="*", label="Observed flux") 
        ax.plot(f_bcens, trans(hist_cube_true[indices][:]), 'r-', marker="*", label="True flux")
    
        for li, label in enumerate(labels):
            xsize = ax.get_xlim()[1] - ax.get_xlim()[0]
            ysize = ax.get_ylim()[1] - ax.get_ylim()[0]
            plt.text(ax.get_xlim()[0] + xsize*0.05, 
                     ax.get_ylim()[0] + ysize*0.10 + li*0.15*ysize , 
                     label, fontsize=10.0)

        ax.set_ylabel('log10(dN/df)')
        ax.set_xlabel('Flux')
        ax.label_outer()

        index = index + 1

    plt.show()

print("Number of sources above max flux in z-bins: ", nabove)

hdus = []
hdus.append(fits.PrimaryHDU(hist_cube))
hdus.append(fits.ImageHDU(data=f_bcens, name="FLUX_OBS_BCENS"))

# This makes it compatible with the old approach
zbins = bin_dict["z_redshift_space_true"]
zbcens = 0.5*(zbins[1:] + zbins[:-1])
hdus.append(fits.ImageHDU(data=zbcens, name="REDSHIFT"))

#for col, bins in iteritems(bin_dict):
#    hdus.append(fits.ImageHDU(data=bins, name=col))

hdul = fits.HDUList(hdus)
hdul.writeto(opts.out)

