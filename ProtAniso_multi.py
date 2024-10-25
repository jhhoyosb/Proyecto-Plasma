# -*- coding: utf-8 -*-

from numpy import linspace, meshgrid
from numpy import array, zeros
from numpy import transpose

# from matplotlib import _cntr as cntr # trace/contourn
import matplotlib.pyplot as plt
from matplotlib.pyplot import contour

from dispeprot import dispeprot_pa
from rootsf import muller_vec, lagrange


def protaniso():
    ## Parameters for instability

    beta_a = .2 # beta plasma relation from temperature and magnetic field
    beta_p = .2
    A_p = 4. # Thermal anisotropy \(T_{i \perp} / T_{i \parallel}\)
    A_a = 4
    gamma_max = 0.5 ## Max expected growth rate _ plot
    n_p = 4.
    n_a = (1 - n_p)/2
    ## Data for initial localization of zeros of D(w,k)=0=Re(D)+i*Im(D)
    # Range for real and imaginary part of frequency
    # X is normalized frequency \(\frac{\omega}{\pm\omega_{ci}}\)
    xr_max = 10
    xr_min = -20
    xi_max = 1
    xi_min = -10
    # Y is normalized wave number \(\frac{ck_z}{\omega_{pi}}\)
    y_max = 4.

    np = 400 # Number of points for each X axis (Re and Im)
    npm1 = np - 1 # Number of intervals defined by NP points
    nir = 2000 # Number of points in Y axis
    nir_interp = 3 # Points without interpolation
    interp_grade = 2 # Polynomial interpolation grade
    epsi = 1.57e-7
    eps_muller = 1e-7

    #nroots = 4 # Number of roots to search

    # Steps for advance in search
    ###dx = (xmax - xmin) / npm1
    ###dy = (ymax - ymin) / npm1

    # Points/mesh to evaluate dispersion function
    xr_ps = linspace(xr_min, xr_max, np)
    xi_ps = linspace(xi_min, xi_max, np)
    Xr_ps, Xi_ps = meshgrid(xr_ps, xi_ps)

    yaxis = linspace(-y_max, y_max, nir)

    #y = y_max # y = yaxis[-1]
    zd = dispeprot_pa(Xr_ps + Xi_ps * 1j, yaxis[-1], beta_p, A_p, n_p, beta_a, A_a, n_a)
    zdr = zd.real
    zdi = zd.imag

    # Plot
    #plt.rc('text', usetex=True)
    fig, ax = plt.subplots(1, 1)
    #ax.hold(True)
    ax.set_xlabel('Frecuencia')
    ax.set_ylabel('Tasa de crecimiento')
    ax.set_title('Contornos de relación de dispersión')
    contour(Xr_ps, Xi_ps, zdr, levels=[0], colors = 'red')
    contour(Xr_ps, Xi_ps, zdi, levels=[0], colors = 'blue')
    # Detect mouse to select root guess
    def onclick(event, roots, figure, axes):
        roots.append([event.xdata, event.ydata])
        if event.button == 3: # right; 1, left
            figure.canvas.mpl_disconnect(cid)
            plt.close(figure)

    roots_guest = []
    on_click = lambda event: onclick(event, roots_guest, fig, ax)
    cid = fig.canvas.mpl_connect('button_press_event', on_click)
    plt.show()
    #plt.draw()
    nroots = len(roots_guest)
    roots = zeros((nir, nroots), dtype=complex)
    roots_yfix = array([guest[0] + guest[1]*1j for guest in roots_guest])
    # First set of roots
    for i in range(nir_interp):
        roots_yfix = muller_vec(roots_yfix, \
                            lambda x: dispeprot_pa(x, yaxis[-(1+i)], beta_p, A_p, n_p, beta_a, A_a, n_a), epsi)
        roots[i,:] = roots_yfix
    roots_yfix = roots[nir_interp-interp_grade-1:nir_interp, :]

    # Rest of roots
    for i in range(nir_interp, nir):
        y_nowi = -(i+1) # Y now index inverse order
        roots_guest = lagrange(array([yaxis[nir - i + interp_grade - k] for k in range(interp_grade + 1)]), roots_yfix, \
                               yaxis[y_nowi])
        roots_guest = muller_vec(roots_guest, \
                            lambda x: dispeprot_pa(x, yaxis[y_nowi], beta_p, A_p, n_p, beta_a, A_a, n_a), eps_muller)
        roots[i, :] = roots_guest
        roots_yfix[0,:] = roots_yfix[1,:]
        roots_yfix[1,:] = roots_yfix[2,:]
        roots_yfix[2,:] = roots_guest
    roots_wr = roots.real
    roots_wi = roots.imag

    fig, ax = plt.subplots(1, 2)
    ax[0].set_xlim([yaxis[0]+1, yaxis[-1]-1])
    ax[0].set_ylim([-2, 2])
    ax[0].set_xlabel('Número de onda')
    ax[0].set_ylabel('Frecuencia$')
    ax[0].plot(yaxis[-1::-1], roots_wr)
    ax[0].set_title('Análisis de dispersión')
    ax[1].set_xlim([0, yaxis[-1]-1])
    ax[1].set_ylim([-1, gamma_max])
    ax[1].set_xlabel('Número de onda')
    ax[1].set_ylabel('Tasa de crecimiento')
    ax[1].set_title('Análisis de tasas de crecimiento')
    ax[1].plot(yaxis[-1::-1], roots_wi)
    # Agregar velocidad de fase respecto a k.
    # Agregar contornos de dispersion
    plt.show()
protaniso()
