
from numpy import linspace, meshgrid, zeros, array
import matplotlib.pyplot as plt
from rootsf import *
from dispeprot import *



class BaseAnisotropy:


    def __init__(self, xr_max, xr_min, xi_max, xi_min, y_max,
                 np, nir, nir_interp, interp_grade, epsi, eps_muller):

        # Data for initial localization of zeros of D(w,k)=0=Re(D)+i*Im(D)
        # Range for real and imaginary part of frequency
        # X is normalized frequency \(\frac{\omega}{\pm\omega_{ci}}\)
        self.xr_max = xr_max
        self.xr_min = xr_min
        self.xi_max = xi_max
        self.xi_min = xi_min
        
        # Y is normalized wave number \(\frac{ck_z}{\omega_{pi}}\)
        self.y_max = y_max

        self.np = np                      # Number of points for each X axis (Re and Im)
        self.nir = nir                    # Number of points in Y axis
        self.nir_interp = nir_interp      # Points without interpolation
        self.interp_grade = interp_grade  # Polynomial interpolation grade
        self.epsi = epsi
        self.eps_muller = eps_muller
        self.yaxis = linspace(-self.y_max, self.y_max, self.nir)


    def create_mesh(self):
        """
        Create the mesh grid for evaluation.
        """

        # Points/mesh to evaluate dispersion function
        xr_ps = linspace(self.xr_min, self.xr_max, self.np)
        xi_ps = linspace(self.xi_min, self.xi_max, self.np)

        return meshgrid(xr_ps, xi_ps)


    def plot_results(self, Xr_ps, Xi_ps, zdr, zdi):
        """
        Common plotting logic for dispersion analysis
        """
        
        # plot dispersion contours
        fig, ax = plt.subplots(1, 1)
        ax.set_xlabel('Frecuencia')
        ax.set_ylabel('Tasa de crecimiento')
        ax.set_title('Contornos de relación de dispersión')

        ax.contour(Xr_ps, Xi_ps, zdr, levels=[0], colors='red')
        ax.contour(Xr_ps, Xi_ps, zdi, levels=[0], colors='blue')

        # detect the mouse clicks for roots selection
        roots_guest: list = []

        def onclick(event):
            roots_guest.append([event.xdata, event.ydata])
            if event.button == 3:  # closing figure with right click
                plt.close(fig)
        
        cid = fig.canvas.mpl_connect('button_press_event', onclick)

        plt.show()

        return roots_guest  # return collected roots for further processing
    

    def plot_roots_analysis(self, yaxis, roots_wr, roots_wi, max_ylim):

        fig, ax = plt.subplots(1, 2)

        ax[0].set_xlim([yaxis[0] + 1, yaxis[-1] - 1])
        ax[0].set_ylim([-2, 2])
        ax[0].set_xlabel('Número de onda')
        ax[0].set_ylabel('Frecuencia')
        ax[0].plot(yaxis[-1::-1], roots_wr)
        ax[0].set_title('Análisis de dispersión')

        ax[1].set_xlim([0, yaxis[-1] - 1])
        
        # Placeholder for gamma_max; this will be set in derived classes
        ax[1].set_ylim([-1, max_ylim])  # Default value; will be overridden. CHANGE THIS
        
        ax[1].set_xlabel('Número de onda')
        ax[1].set_ylabel('Tasa de crecimiento')
        ax[1].set_title('Análisis de tasas de crecimiento')
        
        ax[1].plot(yaxis[-1::-1], roots_wi)
        
        plt.show()
    

    def _find_roots(self, roots_guest, dispeprot_func, **kwargs):
        """
        Common root finding logic. The **kwargs correspond to the arguments that will go on the
        respective dispeprot function
        """

        nroots = len(roots_guest)
        roots = zeros( (self.nir, nroots), dtype=complex )
        roots_yfix = array( [guest[0] + guest[1] * 1j for guest in roots_guest] )

        # first set of roots
        for i in range(self.nir_interp):
            roots_yfix = muller_vec( roots_yfix,
                                     lambda x: dispeprot_func(x, y=self.yaxis[-(1+i)], **kwargs),
                                     self.epsi )
            roots[i, :] = roots_yfix
        
        roots_yfix = roots[ self.nir_interp - self.interp_grade - 1 : self.nir_interp, : ]

        # rest of the roots
        for i in range(self.nir_interp, self.nir):

            y_nowi = - (i + 1)

            roots_guest = lagrange(
                array( [ self.yaxis[self.nir - i + self.interp_grade - k] \
                         for k in range(self.interp_grade + 1) ] ),
                roots_yfix,
                self.yaxis[y_nowi]
            )

            roots_guest = muller_vec( roots_guest, 
                                      lambda x: dispeprot_func(x, y=self.yaxis[y_nowi], **kwargs),
                                      self.eps_muller )
            
            roots[i, :] = roots_guest
            roots_yfix[0, :] = roots_yfix[1, :]
            roots_yfix[1, :] = roots_yfix[2, :]
            roots_yfix[2, :] = roots_guest

        return roots.real, roots.imag


