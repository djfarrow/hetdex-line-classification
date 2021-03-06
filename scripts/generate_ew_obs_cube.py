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
from line_classifier.lfs_ews.equivalent_width import EquivalentWidthAssigner

   
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
    ew = EquivalentWidthAssigner.from_config(config, 'OII_EW')
else:
    ew = EquivalentWidthAssigner.from_config(config, 'LAE_EW')


# Set up equivalent width bins
# Doesn't work well with the fast interpolation methods
#ew_bins_linear = linspace(0.0, 50, 20)
#ew_bins_log = logspace(log10(50 + ew_bins_linear[1] - ew_bins_linear[0]), log10(ew_max), 30)
#ew_bins = concatenate((ew_bins_linear, ew_bins_log))
ew_min = 0.3
ew_bins = logspace(log10(ew_min), log10(ew_max), 50)


ew_bcens = 0.5*(ew_bins[:-1] + ew_bins[1:])
ew_bsize = ew_bins[1:] -ew_bins[0:-1]

nbins = 12

# Dictionary to store the bins for the histogram 
# with keys that are the parameter column name
bin_dict = OrderedDict()

# Redshift range over which OII emitters 
# are actually detectable
if opts.oii:
    # Binsz smaller for OII. ### Want this?
    bin_dict["z_redshift_space_true"] = linspace(0.05, 0.5, nbins)
    #bin_dict["FLUX_OBS"]  = logspace(log10(2e-17), log10(2e-15), nbins)
else:
    bin_dict["z_redshift_space_true"] = linspace(2.18, 3.51, nbins)
    #bin_dict["FLUX_OBS"]  = logspace(log10(2e-17), log10(2e-15), nbins)

# Contruct a tuple with the correct number of
# indices given the number of parameters
# to bin by
shape = []
index_arrays = []
for col, bins in iteritems(bin_dict):
    shape.append(len(bins) - 1)
    index_arrays.append(range(len(bins)))
shape.append(len(ew_bins) - 1)

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

    ew_obs = table["EW_OBS"]
    ew_true = table["EW_TRUE"]
    z = table["z_redshift_space_true"]

    ew_obs_rest = ew_obs/(1.0 + z)
    ew_true_rest = ew_true/(1.0 + z)

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
        indices_in = len(ew_obs_rest)*[True]
        for n, ind in enumerate(indices):
            indices_in = indices_in&(indices_cat[n] == (ind + 1))
            
        tew_obs_rest = ew_obs_rest[indices_in]
        hist_obs = histogram(tew_obs_rest, bins=ew_bins)[0]
 
        hist_cube[indices][:] += hist_obs

        # also include stuff below the range in this number
        nabove[indices] += len(tew_obs_rest[(tew_obs_rest > ew_max)|(tew_obs_rest < ew_min)])

        # Store this for a plot later
        if opts.plot:    
            tew_true_rest = ew_true_rest[indices_in]
            hist_true = histogram(tew_true_rest, bins=ew_bins)[0]
            hist_cube_true[indices] += hist_true


# Loop over doing normalisation
for indices, junk in ndenumerate(zeros(shape[:-1])):   
 
    # Remember to include sources outside of range in normalisation
    norm_obs = (sum(hist_cube[indices][:]) + nabove[indices])*ew_bsize
    frac_under = 1.0*sum(hist_cube[indices][:])/(sum(hist_cube[indices][:]) + nabove[indices])

    hist_cube[indices][:] /= norm_obs
 
    print(frac_under)

    # Avoid NaNs like the plague
    hist_cube[indices][norm_obs < 1e-99] = 0
 


# Plot the EW distributions
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
     
        norm_true = sum(hist_cube_true[indices][:])*ew_bsize
        hist_cube_true[indices][:] /= norm_true
      
        # Base these of the last catalogue used (is ther a better approach??)
        ew_func = ew.ew_function(zmid)  
    
        # Plot of log10 scale
        trans = lambda x : log10(x)
    
        ax = plt.subplot(5, 5, index)
        ax.plot(ew_bcens, trans(hist_cube[indices][:]), 'k-', marker="*", label="Observed EW") 
        ax.plot(ew_bcens, trans(hist_cube_true[indices][:]), 'r-', marker="*", label="True EW")
        ax.plot(ew_bcens, trans(ew_func(ew_bcens)), 'r:', linewidth=3.0, label="Input EW")
    
    
        for li, label in enumerate(labels):
            xsize = ax.get_xlim()[1] - ax.get_xlim()[0]
            ysize = ax.get_ylim()[1] - ax.get_ylim()[0]
            plt.text(ax.get_xlim()[0] + xsize*0.05, 
                     ax.get_ylim()[0] + ysize*0.10 + li*0.15*ysize , 
                     label, fontsize=10.0)

        ax.set_ylabel('log10(dN/dEW)')
        ax.set_xlabel('Equivalent Width')
        ax.label_outer()

        index = index + 1

    plt.show()

print("Number of sources above max EW in z-bins: ", nabove)

hdus = []
hdus.append(fits.PrimaryHDU(hist_cube))
hdus.append(fits.ImageHDU(data=ew_bcens, name="EW_BCENS"))

# This makes it compatible with the old approach
zbins = bin_dict["z_redshift_space_true"]
zbcens = 0.5*(zbins[1:] + zbins[:-1])
hdus.append(fits.ImageHDU(data=zbcens, name="REDSHIFT"))

#for col, bins in iteritems(bin_dict):
#    hdus.append(fits.ImageHDU(data=bins, name=col))

hdul = fits.HDUList(hdus)
hdul.writeto(opts.out)

