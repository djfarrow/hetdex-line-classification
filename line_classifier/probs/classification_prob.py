"""

Carry out a Bayesian classification of a line source based off
off the proposed classifier from Farrow+ (in prep.)

Author: Daniel Farrow 2018 (parts adapted from Andrew Leung's code from Leung+ 2017)


"""
from __future__ import absolute_import

import logging
from numpy import any as nany
from numpy import pi, square, exp, array, power, zeros, ones, isnan, sqrt, mean, append, arange
from scipy.stats import norm

from line_classifier.lfs_ews.luminosity_function import LuminosityFunction
from line_classifier.lfs_ews.equivalent_width import EquivalentWidthAssigner, InterpolatedParameter
from line_classifier.misc.tools import generate_cosmology_from_config

from astropy.table import Table

logging.basicConfig(level=logging.INFO)
_logger = logging.getLogger("lae_prob")

class TooFaintForLimitsException(Exception):
    pass

class NegativeProbException(Exception):
    pass

class UnrecognizedSourceException(Exception):
    pass

def return_delta_volume(wl, lambda_, cosmo, wl_lae=1215.67):
    """
    Return the z*dV/dz at a given wavelength. Scales
    all redshift ranges to the LAE redshift range

    Parameters
    ----------
    wl : array
        the observed wavelengths
    lambda_ : float
        the REST FRAME wavelength of the line to consider
    cosmo : astropy.cosmology:FLRW
        a astropy cosmology object
    wl_lae : float
        the wavelength if Lyman-alpha
        emission (Optional)

    Returns
    -------
    dvdz_delta_z : array
        dV/dz

    """

    zs = wl/lambda_ - 1.0

    dvdz = cosmo.differential_comoving_volume(zs)

    x = cosmo.comoving_distance(zs).to("Mpc")

    # Fixed delta wavelength
    return dvdz*wl_lae/lambda_

def return_lf_n(flux, zs, lfa, cosmo):
    """
 
    Return (L/L*)*dN/dL for a redhsift of
    zs, for an object of a given flux

    Parameters
    ----------
    flux : array
       the fluxes of the sources
    zs : array
       the redshifts of the sources
    lfa : LuminosityFunction
       class to return the LFs
    cosmo : astropy.cosmology
       class for the cosmology
 
    Returns
    -------
    dN/dL : array
       values of dN/dL for the input
       redshifts and luminosity correspoding
       to the input fluxes

    """

    d = cosmo.luminosity_distance(zs).to('cm').value
    ls = flux*(4*pi*d*d)

    lfs = lfa.return_lf_at_z(zs)

    # As we want phi(L)dL = schechter d(L/L*)
    lf_here =  lfs(ls)*(ls/lfa.Lstars_func(zs))

    if lfa.ew_assigner:
        lf_here *= lfa.ew_assigner.classification_correction(zs)

    return lf_here
    

def return_ew_n(ew_obs, zs, ew_func, interp_ew=True):
    """
 
    Return ew*dN/d(EW) for a redhsift of
    zs, for an object of flux flux

    Parameters
    ----------
    ew_obs : array
       the observed equivalent widths
    zs : array
       the redshifts of the sources
    ew_func : EWFunction
       class to assign the equivalent
       widths
    """
    ew_rest = ew_obs/(1. + zs)

    if interp_ew:
        new = ew_func.return_n(zs, ew_rest)*ew_rest   
    else:
        new = (ew_func.ew_function(zs)(ew_rest))*ew_rest   

    return new

def lf_n_interp(fluxes, zs, lfa, finter, cosmo):
    
    n = zeros(len(fluxes)) 

    # For bright stuff just use the LF
    bright = fluxes > finter.maxp
    n[bright] = return_lf_n(fluxes[bright], zs[bright], lfa, cosmo)

    # For faint stuff use interpolated cube
    faint = fluxes <= finter.maxp
    n[faint] = finter.return_new(zs[faint], fluxes[faint])

    return n


def n_additional_line(line_fluxes, line_flux_errors, addl_fluxes, addl_fluxes_error, rel_line_strength, ignore_noise=False):
    """
    Return flux*dN/d(flux) for the flux of an additional
    emission line detected in the spectrum assuming the 
    source is LAE or OII. For LAEs no other lines
    are expected, so it's just this integral from pure
    noise. 

    XXX TODO Shouldn't there be additional scatter for the OIIs
             due to the (possible?) scatter in relative line flux 
             relations?

    Parameters
    ----------
    line_fluxes, line_flux_errors : array or float
        flux(es) and error(s) of the LAE or OII candidate
    addl_fluxes, addl_fluxes_error : array or float
        flux(es) and error(s) the other line
    rel_line_strength : float
        the expected strength of the line relative to OII,
        assuming the source is OII
 
    Returns
    -------
    ndata_lae, ndata_oii : array
        flux*dN/d(flux) for the flux
        of the other emission lines. For 
        LAE and OII
    """ 

    # Error on additional line + error on prediction using relative line strengths
    stddevs = addl_fluxes_error

    if ignore_noise:
        stddevs_with_line = stddevs
    else:
        stddevs_with_line = sqrt(square(stddevs) + square(rel_line_strength*line_flux_errors))

    expected_flux = rel_line_strength*line_fluxes

    ndata_oii = norm.pdf(addl_fluxes, loc=expected_flux,  scale=stddevs_with_line)*abs(addl_fluxes)

    # No line here as this is chance of this flux from just noise
    ndata_lae = norm.pdf(addl_fluxes, loc=0.0,  scale=stddevs)*abs(addl_fluxes)

    # Leave overall probability unchanged if line out of range or not measured
    out_of_range = (stddevs < 1e-30) | (addl_fluxes < -98)

    # Array or floats passed?
    try:
        len(addl_fluxes)
    except TypeError as e:
        if out_of_range:
            ndata_lae = 1.0
            ndata_oii = 1.0
    else:
        ndata_lae[out_of_range] = 1.0 
        ndata_oii[out_of_range] = 1.0

    return ndata_lae, ndata_oii



def source_prob(config, ra, dec, zs, fluxes, flux_errs, ews_obs, ew_err, c_obs, which_color, 
                addl_fluxes, addl_fluxes_error, addl_line_names, flim_file, h=0.67, 
                ignore_noise=False, extended_output=False):
    """
    Return P(LAE|DATA)/P(DATA) and P(DATA|LAE)P(LAE)/(P(DATA|OII)P(OII)) given 

    input information about the sources

    Parameters
    ----------
    config : ConfigParser
        configuration object 
    ra, dec, zs, fluxes, flux_errs, ews_obs, ew_err, c_obs : array
        positions, redshifts (assuming LAE), line fluxes, errors on
        line fluxes, equivalent widths, errors on EW and colours 
        of the sources. (XXX Latter two not used!)
    which_color : str
        which colour is given. Colour not used! 
    addl_fluxes, addl_fluxes_error : 2D array
        each i of addl_fluxes[i, :] should correspond,
        in order, to the fluxes measured for each source
        for the emission line named at position i in
        addl_line_names. To not use set to None (not an
        array of None!)
    addl_line_names : array
        names of emission lines stored in addl_fluxes, in the 
        correct order (see above). To not use set to None 
        (not an array of None!)
    flim_file : str
        file containing the flux limits (here for
        compatibility with Leung+ 2016 style API, not
        used here) 
    h : float
        Hubbles constant/100
    ignore_noise : bool
        ignore the noise on the input
        parameters and assume they are
        perfect
    extended_output : bool
        Return extended output


    Returns
    -------
    prob_lae_given_data :
        probability source is an LAE
    prob_lae_given_data_justlum :
        if Extended = True, probability computed
        just from flux and redshift
    prob_lae_given_data_lum_ew
        if Extended = True, probability computed
        just from flux, redshift and equivalent 
        width 
    prob_lae_given_data_lum_lines:
        if Extended = True, probability computed
        just from flux, redshift and the flux
        in other emission lines
    """

    lae_ew = EquivalentWidthAssigner.from_config(config, 'LAE_EW')
    oii_ew = EquivalentWidthAssigner.from_config(config, "OII_EW")
    lf_lae = LuminosityFunction.from_config(config, "LAE_LF", ew_assigner=lae_ew)
    lf_oii = LuminosityFunction.from_config(config, "OII_LF")
    cosmo = generate_cosmology_from_config(config)

    oii_zlim = config.getfloat("General", "oii_zlim")

    _logger.info("Using Hubbles Constant of {:f}".format(h*100))

    interp_ew = True
    lae_ew_obs = InterpolatedParameter(config.get("InterpolatedEW", "lae_file"), "EW_BCENS")
    oii_ew_obs = InterpolatedParameter(config.get("InterpolatedEW", "oii_file"), "EW_BCENS")

    if ignore_noise:
        lae_ew_obs = lae_ew
        oii_ew_obs = oii_ew
        interp_ew = False

    oii_ew_max = config.getfloat("InterpolatedEW", "oii_ew_max")

    # Cast everything to arrays
    ra = array(ra)
    dec = array(dec)
    zs = array(zs)
    fluxes = array(fluxes)
    ews_obs = array(ews_obs)
    c_obs = array(c_obs)

    wls = (zs + 1.0)*config.getfloat("wavelengths", "LAE")
    zs_oii = wls/config.getfloat("wavelengths", "OII") - 1.0

    # Compute the volume elements 
    dvol_lae = return_delta_volume(wls, config.getfloat("wavelengths", "LAE"), cosmo, wl_lae=config.getfloat("wavelengths", "LAE"))
    dvol_oii = return_delta_volume(wls, config.getfloat("wavelengths", "OII"), cosmo, wl_lae=config.getfloat("wavelengths", "LAE") )

    # Remove source too close to be mistaken for LAEs and therefore removed
    # from catalogue (follows argument used for Leung+ 2017)
    dvol_oii[zs_oii < oii_zlim] = 0.0

    # EW factors
    ew_n_lae = return_ew_n(ews_obs, zs, lae_ew_obs, interp_ew = interp_ew)
    ew_n_oii = return_ew_n(ews_obs, zs_oii, oii_ew_obs, interp_ew = interp_ew)

    # Always LAE according to Leung+ (seems to be true in the sims, very, very rare for OII)
    # Might need to change if EW_ERR grows for OII
    ew_n_oii[ews_obs < 0.0] = 0.0
    ew_n_lae[ews_obs < 0.0] = 1.0

    # Upper limit of the OII EW tabulation - assume LAE
    ew_n_oii[ews_obs > oii_ew_max] = 0.0
    ew_n_lae[ews_obs > oii_ew_max] = 1.0

    # Luminosity function factors
    if ignore_noise:
        lf_n_lae = return_lf_n(fluxes, zs, lf_lae, cosmo)
        lf_n_oii = return_lf_n(fluxes, zs_oii, lf_oii, cosmo)
    else:
        fluxes_lae = InterpolatedParameter(config.get("InterpolatedFlux", "lae_file"), "FLUX_OBS_BCENS")
        fluxes_oii = InterpolatedParameter(config.get("InterpolatedFlux", "oii_file"), "FLUX_OBS_BCENS")

        lf_n_lae = lf_n_interp(fluxes, zs, lf_lae, fluxes_lae, cosmo) 
        lf_n_oii = lf_n_interp(fluxes, zs_oii, lf_oii, fluxes_oii, cosmo)


    # Add additional lines to classification probability (if they're there)
    n_lines_lae = 1.0
    n_lines_oii = 1.0
    if type(addl_line_names) != type(None):
        for line_name, taddl_fluxes, taddl_fluxes_errors in zip(addl_line_names, addl_fluxes[:], addl_fluxes_error[:]):
 
            tn_lines_lae, tn_lines_oii = n_additional_line(fluxes, flux_errs, taddl_fluxes, taddl_fluxes_errors, 
                                                           config.getfloat("RelativeLineStrengths", line_name), 
                                                           ignore_noise=ignore_noise)

            if nany(tn_lines_lae < 0.0) or nany(tn_lines_oii < 0.0):
                dodgy_is = tn_lines_lae < 0.0
                _logger.warning(tn_lines_lae[dodgy_is], fluxes[dodgy_is], taddl_fluxes[dodgy_is], zs_oii[dodgy_is], line_name)
                _logger.warning("The line {:s} results in some negative probabilities".format(line_name))

            # Not an LAE or an OII?
            neither = (tn_lines_lae + tn_lines_oii) < 1e-30
            if nany(neither):
                _logger.warning("Source is neither OII or LAE based off of other emission lines")

            n_lines_lae *= tn_lines_lae
            n_lines_oii *= tn_lines_oii

    # Compute expected number based off of L and z
    nlae = lf_n_lae*dvol_lae
    noii = lf_n_oii*dvol_oii
    prob_lae_given_data_justlum = nlae/(nlae + noii)

    # Include EW information but not emission lines
    nlae_ew = nlae*ew_n_lae
    noii_ew = noii*ew_n_oii

    prob_lae_given_data_lum_ew = nlae_ew/(nlae_ew + noii_ew)

    # Include emission lines but not EW
    nlae_lines = nlae*n_lines_lae
    noii_lines = noii*n_lines_oii
    prob_lae_given_data_lum_lines = nlae_lines/(nlae_lines + noii_lines)

    # Include everything (default return)
    nlae = nlae*n_lines_lae*ew_n_lae
    noii = noii*n_lines_oii*ew_n_oii
    prob_lae_given_data = nlae/(nlae + noii)

    if not extended_output:
        return prob_lae_given_data
    else:
        #return prob_lae_given_data, prob_lae_given_data_justlum, prob_lae_given_data_lum_ew, prob_lae_given_data_lum_lines, ew_n_lae, ew_n_oii 
        return prob_lae_given_data, prob_lae_given_data_justlum, prob_lae_given_data_lum_ew, prob_lae_given_data_lum_lines
