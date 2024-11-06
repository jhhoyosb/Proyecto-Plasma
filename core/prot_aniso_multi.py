
from core.base_anisotropy import BaseAnisotropy
from dispeprot import *



class ProtAnisoMulti(BaseAnisotropy):


    def __init__(self, beta_a, A_a, n_a,
                 beta_p, A_p, n_p,
                 gamma_max, **kwargs):

        # initializing base anisotropy class with common parameters
        super().__init__(**kwargs)

        self.beta_a = beta_a
        self.A_a = A_a
        self.n_a = n_a
        self.beta_p = beta_p
        self.A_p = A_p
        self.n_p = n_p
        self.gamma_max = gamma_max

    
    def run_analysis(self):

        zd = dispeprot_pa(self.Xr_ps + self.Xi_ps * 1j, self.yaxis[-1],
                          self.beta_p, self.A_p, self.n_p,
                          self.beta_a, self.A_a, self.n_a)
        zdr = zd.real
        zdi = zd.imag

        self.plot_results(zdr, zdi)
