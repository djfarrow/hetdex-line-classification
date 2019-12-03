"""

Test implementation of Farrow+ 2018 style probabilities

"""
import pytest
import os
from line_classifier.probs.classification_prob import *
from line_classifier.misc.tools import generate_cosmology_from_config, read_flim_file

def return_lf(config, request, **kwargs):
    """ Tool to allow looping over LF functions """
    if request == "oii":
        return LuminosityFunction.from_config(config, "OII_LF", **kwargs)
    elif request == "lae":
        return LuminosityFunction.from_config(config, "LAE_LF", **kwargs)
    else:
        raise Exception("Unrecognized LF type")

@pytest.fixture
def cosmo(config):
    return generate_cosmology_from_config(config)

@pytest.fixture
def interp_lae_ew(datadir, config):

    fn = datadir.join("config").join(config.get("InterpolatedEW", "lae_file")).strpath

    return InterpolatedParameter(fn, "EW_BCENS")

"""
All these tests just comparisons to results of earlier runs
of the code. So only testing for changes!
"""

@pytest.mark.parametrize("wl, lambda_, expt", [ 
                                               (3900.0, 1215.67, 12789262568.3),
                                               (4800.0, 1215.67, 12609429033.0),
                                               (5200.0, 1215.56, 12320454617.4)
                                              ])

def test_return_delta_volume(cosmo, wl, lambda_, expt):

    dvdz_over_lambda = return_delta_volume(wl, lambda_, cosmo)

    assert dvdz_over_lambda.value == pytest.approx(expt, rel=0.0001)


@pytest.mark.parametrize("type_, flux, zs, expt", [ 
                                                   ("lae", 3e-18, 3.1, 0.01509791),
                                                   ("lae", 1e-17, 1.9, 0.01483313),
                                                   ("oii", 6e-18, 0.24, 0.0372029),
                                                   ("oii", 2e-17, 0.11, 0.03812169)
                                                  ])
def test_return_lf_n(config, cosmo, type_, flux, zs, expt):

    lfa = return_lf(config, type_)
    n = return_lf_n(flux, zs, lfa, cosmo)


    assert n == pytest.approx(expt, rel=0.0001)



@pytest.mark.parametrize("ew_obs, zs, expt", [
                                              (20.0, 2.6, 0.07114728742635372),
                                              (100.0, 3.1, 0.20837608796897072), 
                                              (5.6, 1.9, 0.047429114143167715)
                                             ])

def test_return_ew_n(interp_lae_ew, ew_obs, zs, expt):

    n = return_ew_n(ew_obs, zs, interp_lae_ew)

    assert n == pytest.approx(expt, rel=0.0001)




@pytest.mark.parametrize("name, line_fluxes, addl_fluxes, zoiis, epx_oii, epx_lae", [("NeIII", array([9e-17]), array([5e-17]), array([0.00450909]), 0.02740571, 6.23422985e-06),
                                                                                     ("H_beta", array([9e-17]), array([6e-17]), array([0.01755467]), 1.40145967, 1.044719866e-13),
                                                                                     ("OIII5007", array([9e-17]), array([1e-17]), array([0.01755467]), 6.56691174e-203, 0.21875186)
                                                                                    ])


def test_n_additional_line(config, flim_file, name, line_fluxes, addl_fluxes, zoiis, epx_oii, epx_lae):

    flims = read_flim_file(flim_file)
    wavelength = config.getfloat("wavelengths", name)*(1.0 + zoiis) 
    errors = array(0.2*flims(wavelength/config.getfloat("wavelengths", "LAE") - 1.0)) 

    n_lines_lae, n_lines_oii = n_additional_line(line_fluxes, array([0.0]), addl_fluxes, errors, 
                                                 config.getfloat("RelativeLineStrengths", name))

    assert n_lines_oii == pytest.approx(epx_oii, rel=0.001)
    assert n_lines_lae == pytest.approx(epx_lae, rel=0.001)


@pytest.mark.parametrize("z, fluxes, ew_obs, addl_fluxes, addl_names, e_prob_lae",
                                                                   [
                                                                    (1.9, 9e-17, 40, None, None, 1.0),
                                                                    (2.48, 9e-17, 40, None, None, 0.90076314314944406),
                                                                    (3.18, 9e-17, 40, None, None, 0.151037909544365),
                                                                    (2.5, 9e-17, 40, [[0.23*9e-17]], ["NeIII"], 0.33647488),
                                                                    (2.6, 9e-17, 40, [[0.558*9e-17]], ["H_beta"], 0.77731764),
                                                                    (2.7, 9e-17, 40, [[7e-17], [9e-17*4.752/1.791]], ["OIII4959", "OIII5007"], 0.63853473)
                                                                   ])

def test_source_prob(tmpdir_with_config, config, z, fluxes, ew_obs, flim_file, addl_fluxes, addl_names, e_prob_lae):
    """
    Test source probability
    """

    # Folder with config files
    os.chdir(tmpdir_with_config)

    if type(addl_names) != type(None):

        # Get errors on additional emission lines from flux limit file
        flims = read_flim_file(flim_file)
        errors = []
        for x in addl_names:
            zoii = (1.0 + z)*config.getfloat("wavelengths","LAE")/config.getfloat("wavelengths", "OII") - 1.0
            lambda_ = config.getfloat("wavelengths", x)*(1.0 + zoii)
            errors.append(0.2*flims(lambda_/config.getfloat("wavelengths","LAE") - 1.0))
        errors = array(errors)

        prob_lae_given_data = source_prob(config, [100.0], [0.0], [z], array([fluxes]), array([0.0]), [ew_obs], [0.0], None, None, 
                                          array(addl_fluxes), errors, array(addl_names), flim_file, ignore_noise = True)
    else:
        prob_lae_given_data = source_prob(config, [100.0], [0.0], [z], array([fluxes]), array([0.0]), [ew_obs], [0.0], None, None, None, None, None,
                                          flim_file, ignore_noise = True)


    if e_prob_lae > 0.0:
        assert abs(e_prob_lae - prob_lae_given_data) < 0.05*e_prob_lae
    else:
        assert prob_lae_given_data < 1e-40



