from scipy.special import erfc
from numpy import sqrt, pi, exp


def z_plasma(z0):
    """
    z0: numpy array like of complex points to evaluate in zeta plasma function
    """
    return 1j * sqrt(pi) * exp(-z0**2) * erfc(-1j * z0)
