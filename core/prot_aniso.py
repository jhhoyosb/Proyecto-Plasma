
from core.base_anisotropy import BaseAnisotropy
from dispeprot import *  # temporary import syntax



class ProtAniso(BaseAnisotropy):


    def __init__(self, beta, A, **kwargs):

        # initializing base class with common parameters
        super().__init__(**kwargs)

        self.beta = beta
        self.A = A


    def run_analysis(self):

        zd = dispeprot(self.Xr_ps + self.Xi_ps * 1j,
                       self.beta, self.yaxis[-1], self.A)
        zdr = zd.real
        zdi = zd.imag
        
        self.plot_results(zdr, zdi)

