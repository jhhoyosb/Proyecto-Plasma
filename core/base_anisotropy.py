
from numpy import linspace, meshgrid
import matplotlib.pyplot as plt



class BaseAnisotropy:


    def __init__(self, xr_max, xr_min, xi_max, xi_min, y_max,
                 np, nir, nir_interp, interp_grade, epsi, epsi_muller):

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
        self.eps_muller = epsi_muller

        self._create_mesh()


    def _create_mesh(self):

        # Points/mesh to evaluate dispersion function
        xr_ps = linspace(self.xr_min, self.xr_max, self.np)
        xi_ps = linspace(self.xi_min, self.xi_max, self.np)
        self.Xr_ps, self.Xi_ps = meshgrid(xr_ps, xi_ps)

        self.yaxis = linspace(-self.y_max, self.y_max, self.nir)


    def plot_results(self, zdr, zdi):

        fig, ax = plt.subplots(1, 1)
        
        ax.set_xlabel('Frecuencia')
        ax.set_ylabel('Tasa de crecimiento')
        ax.set_title('Contornos de relación de dispersión')

        plt.contour(self.Xr_ps, self.Xi_ps, zdr, levels=[0], colors='red')
        plt.contour(self.Xr_ps, self.Xi_ps, zdi, levels=[0], colors='blue')

        roots_guest = []
        on_click = lambda event: _onclick(event, roots_guest, fig, ax)
        cid = fig.canvas.mpl_connect('button_press_event', on_click)
        plt.show()

        # Detect mouse to select root guess
        def _onclick(event, roots, figure, axes):
            roots.append([event.xdata, event.ydata])
            if event.button == 3:  # right; 1, left
                figure.canvas.mpl_disconnect(cid)
                plt.close(figure)

        roots_guest = []
        on_click = lambda event: _onclick(event, roots_guest, fig, ax)
        cid = fig.canvas.mpl_connect('button_press_event', on_click)
        
        plt.show()


