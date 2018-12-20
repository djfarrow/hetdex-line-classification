""" Test the equivalent width module """

import pytest

# Compare to an older version of the code (not a guarantee is 100% correct,
# but at least it gives the same results as an older version!). Change
# this if you change the EW in universe.cfg in the config/data dir
@pytest.mark.parametrize( "z, EW, exp", [(0.2, 10, 0.036783678545586909),     
                                         (0.1, 14.2, 0.021185431187433764),
                                         (0.4, 1.25, 0.054294599697860282)])
def test_ew_func(oii_ew_assigner, z, EW, exp):
    """ Test the EW function amplitudes """

    val = oii_ew_assigner.ew_function(z)(EW)
    assert val == pytest.approx(exp)

#
# Test compared to older code version. z=2.1 and z=3.1
# results also published in Leung+ 2017. Weirdly
# I get val=1.47 not 1.49 for z=2.1
@pytest.mark.parametrize("z, exp", [(2.1, 1.471373449240722),
                                    (2.6, 1.3019244496280602),
                                    (3.1, 1.2218485421859402),
                                    (3.45, 1.1870657093033579)]) 
def test_classification_correction(lae_ew_assigner, z, exp):
    """
    Test the classification correction, the correction
    to get total number of sources from measurement that
    assumes everything with EW<20A doesn't belong to 
    the sample.
    """
    val = lae_ew_assigner.classification_correction(z)
    assert val == pytest.approx(exp)

