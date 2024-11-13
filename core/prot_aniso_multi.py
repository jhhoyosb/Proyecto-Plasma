
from core.base_anisotropy import BaseAnisotropy
from dispeprot import dispeprot_pa



class ProtAnisoMulti(BaseAnisotropy):
    
    """
    Must keep in mind that the args for dispeprot_pa function are:
        beta_p, A_p, n_p, beta_a, A_a, n_a
    """

    def __init__(self, beta_a, A_a,
                 beta_p, A_p, n_p,
                 gamma_max, **kwargs):

        # initializing base anisotropy class with common parameters
        super().__init__(**kwargs)

        self.beta_a = beta_a
        self.A_a = A_a
        self.beta_p = beta_p
        self.A_p = A_p
        self.n_p = n_p
        self.n_a = (1 - self.n_p) / 2
        self.gamma_max = gamma_max

    
    def run_analysis(self):

        # creating mesh
        Xr_ps, Xi_ps = self.create_mesh()

        zd = dispeprot_pa(Xr_ps + Xi_ps * 1j, self.yaxis[-1],
                          self.beta_p, self.A_p, self.n_p,
                          self.beta_a, self.A_a, self.n_a)
        zdr = zd.real
        zdi = zd.imag

        # calling common plot method and get roots guest from user input
        roots_guest = self.plot_results(Xr_ps, Xi_ps, zdr, zdi)

        # perform root finding using common method
        roots_wr, roots_wi = self._find_roots(
            roots_guest, dispeprot_pa,
            beta_p=self.beta_p, A_p=self.A_p, n_p=self.n_p,
            beta_a=self.beta_a, A_a=self.A_a, n_a=self.n_a )

        self.plot_roots_analysis(self.yaxis, roots_wr, roots_wi, self.gamma_max)

