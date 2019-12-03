"""

Add LAE/OII probability to to an existing catalogue

"""

import argparse
from configparser import RawConfigParser
from numpy import array
from astropy.table import Table
from line_classifier.probs.classification_prob_leung import source_prob as source_prob_leung
from line_classifier.probs.classification_prob import source_prob 


parser = argparse.ArgumentParser(description="Add P(LAE) to a catalogue")
parser.add_argument("fconfig", help="Config file")
parser.add_argument("filename", help="Name of file to add P(LAE) to")
parser.add_argument("ofilename", help="Name of file to output to")
#parser.add_argument("--h0", help="H(z=0)/100", default=1.0)
opts = parser.parse_args()

# Read in config
config = RawConfigParser()
config.read(opts.fconfig)

# Other emission lines to include, set to None to not use this
addl_el_names = ["NeIII", "H_beta", "OIII4959", "OIII5007"]
#addl_el_names = None

table = Table.read(opts.filename)

# Prepare suitable lists of the flux of other emission lines
# and their names
if type(addl_el_names) != type(None):

    addl_fluxes = array([table[x] for x in addl_el_names])
    addl_fluxes_error = array([table[x + "_ERROR"] for x in addl_el_names])


else:
    addl_fluxes = None
    addl_fluxes_error = None


#
#
#

table["PLAE_noerr"], table["PLAE_lum_noerr"], table["PLAE_ew_noerr"], table["PLAE_lines_noerr"] = source_prob(config, table["ra"], table["dec"], table["z"], table["flux"], table["FLUX_ERR"], 
                                                                                                              table["EW_TRUE"], table["EW_ERR"], [None]*len(table), None, addl_fluxes,
                                                                                                              addl_fluxes_error, addl_el_names, "../common_configuration_stuff/Line_flux_limit_5_sigma_baseline.dat", 
                                                                                                              extended_output=True, ignore_noise=True)


table["PLAE_ignore_errors"], table["PLAE_lum_ignore"], table["PLAE_ew_ignore"], table["PLAE_lines_ignore"] = source_prob(config, table["ra"], table["dec"], table["z"], table["FLUX_OBS"], table["FLUX_ERR"], 
                                                                                                                         table["EW_OBS"], table["EW_ERR"], [None]*len(table), None, addl_fluxes,
                                                                                                                         addl_fluxes_error, addl_el_names, "../common_configuration_stuff/Line_flux_limit_5_sigma_baseline.dat", 
                                                                                                                         extended_output=True, ignore_noise=True)


#table["PLAE_no_flux_err"], table["PLAE_lum_no_flux_err"], table["PLAE_ew_no_flux_err"], table["PLAE_lines_no_flux_err"] = source_prob(config, table["ra"], table["dec"], table["z"], table["flux"], table["FLUX_ERR"], 
#                                                                                                                                      table["EW_OBS"], table["EW_ERR"], [None]*len(table), None, addl_fluxes,
#                                                                                                                                      addl_fluxes_error, addl_el_names, "../common_configuration_stuff/Line_flux_limit_5_sigma_baseline.dat", 
#                                                                                                                                      extended_output=True, ignore_noise=True) 

#throw_away, table["PLAE_leung_noerr"] = source_prob_leung(config, table["ra"], table["dec"], table["z"], table["flux"], table["FLUX_ERR"],
#
#                                                          table["EW_TRUE"], table["EW_ERR"], [None]*len(table), None, addl_fluxes,
#                                                          addl_fluxes_error, addl_el_names, "../common_configuration_stuff/Line_flux_limit_5_sigma_baseline.dat")


#throw_away, table["PLAE_leung_biggerrange"] = source_prob_leung(config, table["ra"], table["dec"], table["z"], table["FLUX_OBS"], table["FLUX_ERR"],
#                                                                table["EW_OBS"], table["EW_ERR"], [None]*len(table), None, addl_fluxes,
#                                                                addl_fluxes_error, addl_el_names, "../common_configuration_stuff/Line_flux_limit_5_sigma_baseline.dat")



table.write(opts.ofilename)
