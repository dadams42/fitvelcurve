import pytest
import numpy as np
import astropy
from fitvelcurve.gal_class import Galaxy

def test_ptform():
    '''
    Check the output of ptform() on a few specific cases
    '''

    # Create an instance of the Galaxy() class
    r_bins = np.linspace(1, 31, 5)
    gal = Galaxy(r_bins)

    # Check both inputs equal to 0 (should return min values of each)
    x = np.array([0, 0])
    params = gal.ptform(x)
    assert params == pytest.approx(np.array([6, 1]))

    # Check both inputs equal to 1 (should return max values for each)
    x = np.array([1, 1])
    params = gal.ptform(x)
    assert params == pytest.approx(np.array([13, 50]))

    # Check a random point in the middle
    x = np.random.uniform(size = 2)
    params = gal.ptform(x)
    # Do the transformation by hand
    true_params = np.zeros(2)
    true_params[0] = (13 - 6) * x[0] + 6
    true_params[1] = (50 - 1) * x[1] + 1
    # Check that they agree
    assert params == pytest.approx(true_params)

def test_rho_nfw():
    '''
    Check the NFW profile in the Galaxy() class vs. the built-in astropy version
    '''

    # Set mvir and cvir
    mvir = 1e13
    cvir = 20
    # Create array of radial points
    r_bins = np.linspace(1, 31, 5)

    # Initiate a Galaxy() object
    test_gal = Galaxy(r_bins, mvir = mvir, cvir = cvir)
    # Compute the NFW density and mass profiles from Galaxy()
    rho_nfw = test_gal.rho_nfw()

    # Compute the same thing with astropy
    astropy_nfw_model = astropy.modeling.physical_models.NFW(mass = mvir, 
                                                             concentration = cvir,
                                                             massfactor = ('critical', 102))
    astropy_rho_nfw = astropy_nfw_model.evaluate(r_bins, mass = mvir, 
                                                concentration = cvir,
                                                redshift = 0.0)

    # Check that these are equal within tolerance
    assert rho_nfw == pytest.approx(astropy_rho_nfw, rel = 1e-2)