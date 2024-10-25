from z_plasma import z_plasma
from numpy import sqrt


def dispeprot(x, beta, y, A):
    # Eq. 227 from Davidson
    # Multiply by -\omega_{pi}^{-2}
    # \(x=\frac{\omega}{\pm\omega_{ci}}\)
    # \(y = \frac{ck_z}{\omega_{pi}}\)
    # \( A = \frac{T_{i \perp}}{T_{i\parallel}}\)
    # \( vth = \frac{Vth \omega_{pi}}{c\omega_{ci}}\)
    vth = sqrt(beta)
    psi = (x - 1.) / (y * vth)
    drst = y**2 + x - (psi * A + 1./(vth * y)) * z_plasma(psi) - A + 1.
    return drst

def dispeprot_pa(x, y, beta_p, A_p, n_p, beta_a, A_a, n_a):
    vth_p = sqrt(beta_p)
    vth_a = sqrt(beta_a)
    psi_p = (x - 1.) / (y * vth_p)
    psi_a = (x - .5) / (y * vth_a)
    dis_p = x - (psi_p * A_p + 1 / (y * vth_p)) * z_plasma(psi_p) - A_p + 1
    dis_a = -2*x - (psi_a*A_a + 1/(y*vth_a))*z_plasma(psi_a) - A_a + 1
    drst = y**2 + dis_p + dis_a * n_a / n_p
    return drst

def dispeprot_multi(x, y, beta, A, n, m, q):
    vth = sqrt(beta)
    psi = (x - q/m)/(y*vth)
    dis = (n*(m*x-(psi*A+1/(y*vth))*z_plasma(psi)-A+1))/n[0]
    drst = y**2 - sum(dis)
    return drst
