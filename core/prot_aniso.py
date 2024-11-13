
from core.base_anisotropy import BaseAnisotropy
from dispeprot import dispeprot  # temporary import syntax



class ProtAniso(BaseAnisotropy):
    
    """
    Must keep in mind that the args for the dispeprot function are:
        beta, A
    """

    def __init__(self, beta, A, **kwargs):

        # initializing base class with common parameters
        super().__init__(**kwargs)

        self.beta = beta
        self.A = A


    def run_analysis(self):

        # creating the mesh
        Xr_ps, Xi_ps = self.create_mesh()

        zd = dispeprot( Xr_ps + Xi_ps * 1j,
                        beta=self.beta, y=self.yaxis[-1], A=self.A )
        zdr = zd.real
        zdi = zd.imag

        # calling common plot method and get roots guest from user input
        roots_guest = self.plot_results(Xr_ps, Xi_ps, zdr, zdi)

        # perform root finding using common method
        roots_wr, roots_wi = self._find_roots(
            roots_guest, dispeprot, beta=self.beta, A=self.A )

        self.plot_roots_analysis(self.yaxis, roots_wr, roots_wi, 0.3)

