"""

Test the luminosity function module.

AUTHORS: Daniel Farrow 2018 (MPE)

"""

from __future__ import absolute_import, print_function

import pytest 
from mpmath import gammainc as gammainc_mpmath
from line_classifier.lfs_ews.luminosity_function import LuminosityFunction, gamma_integral_limits

def return_lf(config, request, **kwargs):
    """ Fixture to allow looping over LF functions """
    if request == "oii": 
        return LuminosityFunction.from_config(config, "OII_LF", **kwargs)
    elif request == "lae":
        return LuminosityFunction.from_config(config, "LAE_LF", **kwargs)
    else:
        raise Exception("Unrecognized LF type")



# Compare to an older version of the code (not a guarantee is 100% correct,
# but at least it gives the same results as an older version!). Change
# this if you change the LF in universe.cfg in the config/data dir
@pytest.mark.parametrize("stype, z, lum, exp", [("lae", 2.156, 1e41, 0.34907584),
                                               ("lae", 2.156, 2e42, 0.00099136),
                                               ("lae", 3.123, 5e40, 1.58085127),
                                               ("oii", 0.213, 5e40, 0.01300613),
                                               ("oii", 0.123, 1e39, 1.84715227)]) 
def test_luminosity_function(config, stype, z, lum, exp):
    """
    Test the luminosity function returns the
    expected values. Compare to an older
    implementation of the LF.

    """
    lf = return_lf(config, stype)
    val = lf.return_lf_at_z(z)(lum)

    # Approximately the same as the expected value
    assert val == pytest.approx(exp)



@pytest.mark.parametrize("s", [-0.9, 0.4, -0.5, -0.65, 2])
@pytest.mark.parametrize("low, high", [(0.01, 1.0), (1.0, 4.0), (0.56, 3.54), (0.1, 1000)])
def test_gamma_integrator(s, low, high):
    """
    Test the gamma function integrator 
    against mpmath
    """
    mpmath = gammainc_mpmath(s, low, high)
    here = gamma_integral_limits([s], [low], [high])

    assert abs(mpmath - here) < 0.001*here


