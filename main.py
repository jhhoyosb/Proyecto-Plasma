# here goes the implementation with typer
import typer

from core.prot_aniso import *
from core.prot_aniso_multi import *


# initializing app
app = typer.Typer()


@app.command()
def run_protaniso( beta: float = 0.2,
                   A: float = 0.5,
                   xr_max: float = 10,
                   xr_min: float = -20,
                   xi_max: float = 1,
                   xi_min: float = -5,
                   y_max: float = 4.,
                   np: int = 400,
                   nir: int = 2000,
                   nir_interp: int = 3,
                   interp_grade: int = 2,
                   epsi: float = 1.57e-7,
                   eps_muller: float = 1e-7 ):

    protaniso_simulator = ProtAniso(
        beta = beta,
        A = A,
        xr_max = xr_max,
        xr_min = xr_min,
        xi_max = xi_max,
        xi_min = xi_min,
        y_max = y_max,
        np = np,
        nir = nir,
        nir_interp = nir_interp,
        interp_grade = interp_grade,
        epsi = epsi,
        eps_muller = eps_muller
    )

    protaniso_simulator.run_analysis()


@app.command()
def run_protaniso_multi( beta_a: float = 0.2,
                         beta_p: float = 0.2,
                         A_a: float = 4.0,
                         A_p: float = 4.0,
                         n_p: float = 4.0,
                         gamma_max: float = 0.5,
                         xr_max: float = 10,
                         xr_min: float = -20,
                         xi_max: float = 1.0,
                         xi_min: float = -10.0,
                         y_max: float = 4.0,
                         np: int = 400,
                         nir: int = 2000,
                         nir_interp: int = 3,
                         interp_grade: int = 2,
                         epsi: float = 1.57e-7,
                         eps_muller: float = 1e-7 ):

    protaniso_multi_simulator = ProtAnisoMulti(
        beta_a = beta_a,
        beta_p = beta_p,
        A_a = A_a,
        A_p = A_p,
        n_p = n_p,
        gamma_max = gamma_max,
        xr_max = xr_max,
        xr_min = xr_min,
        xi_max = xi_max,
        xi_min = xi_min,
        y_max = y_max,
        np = np,
        nir = nir,
        nir_interp = nir_interp,
        interp_grade = interp_grade,
        epsi = epsi,
        eps_muller = eps_muller
    )

    protaniso_multi_simulator.run_analysis()


if __name__ == '__main__':
    app()
