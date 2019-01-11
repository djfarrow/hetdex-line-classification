"""

Test the implementation of Leung+ 2017 style probabilities

"""

import pytest
from numpy import pi, square, linspace, trapz, array
from line_classifier.probs.classification_prob_leung import (LFIntegrator, prob_additional_line, source_prob,
                                                             luminosity_likelihoods, ew_prob)
from line_classifier.lfs_ews.luminosity_function import LuminosityFunction
from line_classifier.misc.tools import generate_cosmology_from_config, read_flim_file

def return_lf(config, request, **kwargs):
    """ Tool to allow looping over LF functions """
    if request == "oii":
        return LuminosityFunction.from_config(config, "OII_LF", **kwargs)
    elif request == "lae":
        return LuminosityFunction.from_config(config, "LAE_LF", **kwargs)
    else:
        raise Exception("Unrecognized LF type")

@pytest.mark.parametrize("case, z", [("lae", 2.5), ("oii", 0.3)])
def test_lf_integrator(config, case, z, flim_file):
    """
    Compare results of the lf_integrator class to
    just doing Trapezium rule on the luminosity
    function

    """

    if case == "lae":
        lambda_ = 1215.67
    elif case == "oii":
        lambda_ = 3727.3
    else:
        assert False 

    lf = return_lf(config, case)
    lf_integrator = LFIntegrator(lf, flim_file)

    lmin = 1e-16*4.0*pi*square(lf.cosmo.luminosity_distance(z).to('cm').value)
    lmax = lf.return_lmax_at_z(z) 

    # huge number of bins to compare to analytic integral 
    bins = linspace(lmin, lmax, 400000) 
    lfmod_unnorm = lf.return_lf_at_z(z, normed=False) 
    n_per_mpc_ref = trapz(lfmod_unnorm(bins), x=bins/lf.Lstars_func(z))

    # compare to lf integrator integral
    cosmo = generate_cosmology_from_config(config)
    integrals, norms = lf_integrator.integrate_lf_limits(config, [100.0], [0.0], [z], [lmin], [lmax], 
                                                         lambda_, cosmo)

    assert abs(integrals[0] - n_per_mpc_ref) < 0.01*n_per_mpc_ref


#
#
# The stuff below is all tested against Andrew Leung's code,
# from Leung+ 2017
#

# 1% of additional line prob from Leung code as results should be near identical
@pytest.mark.parametrize("name, line_fluxes, addl_fluxes, zoiis, epx_oii, epx_lae", [("NeIII", array([9e-17]), array([5e-17]), array([0.00450909]), 0.0029627, 7.945545488e-07),
                                                                                     ("H_beta", array([9e-17]), array([6e-17]), array([0.01755467]), 0.1420576, 1.044719866e-13),
                                                                                     ("OIII5007", array([9e-17]), array([1e-17]), array([0.01755467]), 1.210356811e-203, 0.02188774)
                                                                                    ])

def test_prob_additional_line(config, name, line_fluxes, addl_fluxes, zoiis, flim_file, epx_oii, epx_lae):

    flims = read_flim_file(flim_file)
    wavelength = config.getfloat("wavelengths", name)*(1.0 + zoiis) 
    errors = array(0.2*flims(wavelength/config.getfloat("wavelengths", "LAE") - 1.0)) 

    plae, poii = prob_additional_line(name, line_fluxes, array([0.0]), addl_fluxes, errors,
                                      config.getfloat("RelativeLineStrengths", name))

    assert plae[0] == pytest.approx(epx_lae, rel=1e-2, abs=0.0)
    assert poii[0] == pytest.approx(epx_oii, rel=1e-2, abs=0.0)


"""
Had to Change Andrew Leung's O_m and O_lambda to match my cosmology. Also
had to rescale Andrew Leung's wavelength versus sensitivity tabulation to correct for
using 1216.0 for the z->wavelength conversion.
"""
@pytest.mark.parametrize("flux, case, z, expected", [ (9e-17, "oii", 0.363298, 0.041268), 
                                                      (9e-17, "lae", 1.9, 0.147877),
                                                      (9e-17, "lae", 3.18, 0.042753), 
                                                      (1.2e-16, "oii", 0.363262, 0.037168),
                                                      (1.2e-16, "lae", 1.9, 0.100235),
                                                      (1.2e-16, "lae", 3.18, 0.022630)
                                                    ])
def test_luminosity_likelihoods(config, flux, case, z, expected, flim_file):
    """
    Compare luminoisty likelihoods to results from Andrew Leungs code
    """
    lf = return_lf(config, case)
    if case == "lae":
        lambda_ = 1215.668
    elif case == "oii":
        lambda_ = 3727.45
    else:
        assert False

    cosmo = generate_cosmology_from_config(config)
    pl, n = luminosity_likelihoods(config, [300.0], [0.0], [z], [flux], lf, lambda_, flim_file, cosmo)


    assert pl > 0.0
    assert n > 0.0
    # Agree with Andrew Leung's code to 5%
    assert pl == pytest.approx(expected, rel=0.05)

#
# Test the probability of EW against Leung+ code
#
@pytest.mark.parametrize("mode, ew_obs, z, expt", [("LAE", 40, 1.8999999999999999, 0.023584160100256324),
                                                   ("LAE", 40, 2.48, 0.013929176697349166),
                                                   ("LAE", 40, 3.1799999999999997, 0.0084181930779243519),
                                                   ("OII", 40, 0.13496482581925973, 0.0072311149164594607),
                                                   ("OII", 40, 0.36326234825416837, 0.028735638637520183)]) 
def test_ew_prob(lae_ew_assigner, oii_ew_assigner, mode, ew_obs, z, expt):

   if mode == "LAE":
       ew_func = lae_ew_assigner
   elif mode == "OII":
       ew_func = oii_ew_assigner

   val_lae = ew_prob(array([ew_obs]), array([z]), ew_func)
   print(ew_obs, z, val_lae, expt)
 
   assert val_lae == pytest.approx(expt, rel=1e-3)


# 10% of posterior odds ratio and 5% of prob_lae compared to original Andrew Leung code
# need to be careful here as parameters not valid for any source will raise exception! 
@pytest.mark.parametrize("z, fluxes, ew_obs, addl_fluxes, addl_names, e_ratio, e_prob_lae",
                                                                   [
                                                                    (1.9, 9e-17, 40, None, None, 1e+32, 1.0),
                                                                    (2.48, 9e-17, 40, None, None, 9.0769011810393501, 0.90076314314944406),
                                                                    (3.18, 9e-17, 40, None, None, 0.17790889751426178, 0.151037909544365),
                                                                    (2.08, 9e-17, 40, [[5e-17]], ["NeIII"], 10.917948575339162, 0.91609294219734949),
                                                                    (2.12, 9e-17, 40, [[6e-17]], ["H_beta"], 2.2721726484396545e-09, 2.2721726536024229e-09),
                                                                    (2.08, 9e-17, 40, [[7e-17], [9e-17*4.752/1.791]], ["OIII4959", "OIII5007"], 0.0, 0.0),
                                                                    (2.373678804418128, 4.355691265144264E-17, 34.99668429707373, 
                                                                     [[-1.338393901007366E-18], [2.0183850784203426E-17], [2.0248272025151368E-17]],
                                                                     ["NeIII", "H_beta", "OIII4959"], 2.4124703809980317, 0.70695716347652693)
                                                                    ])
                                                                   
def test_source_prob(config, z, fluxes, ew_obs, flim_file, addl_fluxes, addl_names, e_ratio, e_prob_lae):
    """
    Test source probability
    """


    if type(addl_names) != type(None):

        # Get errors on additional emission lines from flux limit file
        flims = read_flim_file(flim_file)
        errors = []
        for x in addl_names:
            zoii = (1.0 + z)*config.getfloat("wavelengths","LAE")/config.getfloat("wavelengths", "OII") - 1.0
            lambda_ = config.getfloat("wavelengths", x)*(1.0 + zoii)
            errors.append(0.2*flims(lambda_/config.getfloat("wavelengths","LAE") - 1.0))
        errors = array(errors)

        posterior_odds, prob_lae_given_data = source_prob(config, [100.0], [0.0], [z], array([fluxes]), array([0.0]), [ew_obs], [0.0], None, None, 
                                                          array(addl_fluxes), errors, array(addl_names), flim_file)
    else:
        posterior_odds, prob_lae_given_data = source_prob(config, [100.0], [0.0], [z], array([fluxes]), array([0.0]), [ew_obs], [0.0], None, None, None, None, None,
                                                          flim_file)

    if e_prob_lae > 0.0:
        assert abs(e_prob_lae - prob_lae_given_data) < 0.05*e_prob_lae
    else:
        assert prob_lae_given_data < 1e-40

    if e_ratio > 1e-40 :
        assert abs(e_ratio - posterior_odds) < 0.1*e_ratio
    else:
        assert posterior_odds < 1e-40


