#from scipy.special import erfc
#from numpy import sqrt, pi, exp
from plasmapy.dispersion import plasma_dispersion_func


def z_plasma(z0) -> complex:
    """
    Calculate the plasma dispersion function.

    Parameters:
    z0: numpy array like of complex points to evaluate in zeta plasma function.

    Returns:
    complex: The value of the plasma dispersion function.
    """
    #return 1j * sqrt(pi) * exp(-z0**2) * erfc(-1j * z0)
    return plasma_dispersion_func(z0)
