import numpy as np
from numpy import pi, log
from astropy.constants import G, M_sun, kpc
from astropy import units as u
from astropy.cosmology import Planck15 as cosmo
import dynesty
from dynesty import plotting as dyplot
###


class Galaxy():
    """
    Derive a dark matter density profile from a galaxy rotation curve

    Args:
        mvir (float): galaxy virial mass (Msun)
        cvir (float): dark matter concentration parameter
        data_vel (array): galaxy velocity curve data (km/s)
        data_err (array): errors on galaxy velocity curve data (km/s)
        data_r (array): radial positions corresponding to galaxy velocity curve data (kpc)

    Attibutes:
        rhocrit (float): critical density of the universe (Msun/kpc^3)
        fac (int): critical overdensity for defining virial radius
        mvir (float): galaxy virial mass (Msun)
        cvir (float): dark matter concentration parameter
        data_vel (array): galaxy velocity curve data (km/s)
        data_err (array): errors on galaxy velocity curve data (km/s)
        data_r (array): radial positions corresponding to galaxy velocity curve data (kpc)
    """

    def __init__(self, data_r, mvir=1e10, cvir=20, data_vel=[], data_err=[]):

        # Check that data is in the correct format
        assert len(data_r) > 0, "An array of radial values must be provided"
        if len(data_vel) > 0:
            assert data_vel.shape == data_r.shape, "The velocity array must have the same length as the radial data"
        if len(data_err) > 0:
            assert data_err.shape == data_vel.shape, "The velocity error array must have the same length as the velocity array"

        # Define constants
        self.fac = 102
        self.rhocrit = (3*(100*cosmo.h/(kpc*u.kpc))**2/(8*np.pi*G) * (kpc*u.kpc)**3 / (M_sun*u.Msun)).value

        # Data arrays
        self.data_vel = data_vel # Velocity
        self.data_err = data_err # Velocity errors
        self.data_r = data_r # Radii

        # Initialize virial mass and concentration
        self.mvir = mvir # Units: solar mass
        self.cvir = cvir # Unitless

    def delc(self):
        """
        Helper integral for computing the mass enclosed by the NFW profile

        Returns:
            float: value of helper integral
        """
        return (self.fac/3)*self.cvir**3/(np.log(1+self.cvir)-self.cvir/(1+self.cvir))

    def rho_nfw(self):
        """
        NFW density (virial mass in Msun, concentration parameter)
        """
    
        rv = (self.mvir/((4/3) * np.pi * self.fac * self.rhocrit))**(1/3) 
        rs = rv / self.cvir
        rhos = self.rhocrit*self.delc()
    
        return rhos/((self.data_r/rs)*(1.+(self.data_r/rs))**2.)

    def mass_nfw(self):
        """
        Enclosed mass (virial mass in Msun, concentration parameter) in Msun

        Returns:
            array: mass profile at radial position of data_r
        """
    
        rv = (self.mvir/((4/3) * np.pi * self.fac * self.rhocrit))**(1/3)
        rs = rv / self.cvir
        rhos = self.rhocrit*self.delc()
    
        return 4 * np.pi * rhos * rs**3 * (np.log((self.data_r + rs) / rs) + rs/(self.data_r + rs) - 1) 
    
    def velocity(self):
        """
        Calculate the velocity rotation curve (in km/s)

        """
        # Calculate enclosed mass array with NFW profile
        enclosed_mass = self.mass_nfw()

        # Calculate predicted velocity from v^2 = GM/r
        pred_vel = np.sqrt(G.value * enclosed_mass * M_sun.value / self.data_r / kpc.value) / 1000
        return pred_vel
    
    def log_like(self, x):
        """
        Calculate log-likelihood using the data and NFW parameters (c, log10(M_h))

        Args:
            x (array): NFW parameters [log10(mvir), cvir]

        Returns:
            float: log likelihood
        """

        # Read parameters from array
        log_m_h, c = x

        # Redefine the class attributes to match these
        self.mvir = 10**log_m_h
        self.cvir = c

        pred_vel = self.velocity()

        # Calculate log-likelihood
        llh = -1/2 * np.sum(((self.data_vel - pred_vel) / self.data_err) ** 2)

        return llh

    def ptform(self, x):
        """
        Prior transform

        Args:
            x (array): two uniform samples between 0 and 1

        Returns:
            array: NFW parameters [log10(mvir), cvir] 
        """

        # Make sure x is a numpy array
        x = np.array(x)

        # Transform x[0] (log_m_h) to range 6-13 in log(M/M_sun)
        x[0] = (13 - 6) * x[0] + 6
        # Transform x[1] (c) to range 1-50
        x[1] = (50 - 1) * x[1] + 1

        # Return transformed x
        return x
    
    def run_sampler(self):
        """
        Runs sampler

        Returns:
            dynesty sampler object
        """
        sampler = dynesty.DynamicNestedSampler(self.log_like, self.ptform, ndim = 2, walks = 50)
        sampler.run_nested()

        ##plotting results:

        resu = sampler.results 

        print("number of iterations = ", resu.niter)


        # Plot a summary of the run.
        rfig, raxes = dyplot.runplot(resu)

        # Plot traces and 1-D marginalized posteriors.
        tfig, taxes = dyplot.traceplot(resu)

        # Plot the 2-D marginalized posteriors.
        cfig, caxes = dyplot.cornerplot(resu)

        return sampler