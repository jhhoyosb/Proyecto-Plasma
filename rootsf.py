import numpy as np
from numpy import sqrt, array


def muller(xi, f, epsi, tol=1e-10):
    if abs(f(xi)) > tol:
        disp = epsi * abs(xi)
        disp = disp if disp > tol else tol
        xim1 = xi - disp
        xim2 = xi + disp
        error = tol + 1
        while error > tol:
            q = (xi - xim1) / (xim1 - xim2)
            A = q*f(xi) - q*(1+q)*f(xim1) + q**2 * f(xim2)
            B = (2*q + 1)*f(xi) - (1+q)**2*f(xim1) + q**2 * f(xim2)
            C = (1+q) * f(xi)
            Dp = B + sqrt(B**2 - 4*A*C)
            Dm = B - sqrt(B**2 - 4*A*C)
            D = Dp if abs(Dp) >= abs(Dm) else Dm
            x_next = xi - (xi - xim1) * 2*C / D
            xim2 = xim1
            xim1 = xi
            xi = x_next
            error = abs(xi - xim1)
    return xi


def muller_vec(x_points, f, epsi, tol=1e-10):
    return array([muller(x, f, epsi, tol) for x in x_points])


def get_lagrange_basis(x_points, x_eval):
    """
    Computes the Lagrange basis polynomials evaluated at x_eval.

    Args:
    x_points: A numpy array of x-coordinates.
    x_eval: The point at which to evaluate the basis polynomials.

    Returns:
    A numpy array of the evaluated basis polynomials.
    """

    # set up size of points and array to store basis polynomials
    n = len(x_points)
    basis_polynomials = np.zeros(n)

    # main loop to get each basis polynomial
    for i in range(n):
        basis_polynomials[i] = 1.0
        for j in range(n):
            if j != i:
                basis_polynomials[i] *= (x_eval - x_points[j]) / (x_points[i] - x_points[j])

    return basis_polynomials


def lagrange(x_points, y_points, x_eval) -> np.ndarray:
    """
    Computes a 1D array of Lagrange polynomials evaluated at given x_eval single points.

    Args:
    x_points (1D array): numpy 1D array-like of x values to compute each polynomial.
    y_points (2D array): numpy 2D array-like, it stores rows of values for y for a vectorized computing.
    x_eval (float): scalar value of x for evaluation of the Lagrange Polynomial.

    Returns:
    1D array of Lagrange Polynomials evaluated at a given x_eval point.
    """

    lagrange_basis = get_lagrange_basis(x_points, x_eval)
    lagrange_poly = np.dot(y_points.T, lagrange_basis)

    return lagrange_poly


# previous implementation
# def lagrange(x_points, y_points, x_eval):

#     # Generalizar despues...
#     l0 = ((x_eval - x_points[1]) / (x_points[0] - x_points[1])) * \
#         ((x_eval - x_points[2]) / (x_points[0] - x_points[2]))
#     l1 = ((x_eval - x_points[0]) / (x_points[1] - x_points[0])) * \
#         ((x_eval - x_points[2]) / (x_points[1] - x_points[2]))
#     l2 = ((x_eval - x_points[0]) / (x_points[2] - x_points[0])) * \
#         ((x_eval - x_points[1]) / (x_points[2] - x_points[1]))

#     result = y_points[0, :] * l0 + y_points[1, :] * l1 + y_points[2, :] * l2
#     print(x_points.shape, y_points.shape, result.shape, x_eval)

#     all_ls = np.array([l0, l1, l2])

#     #return y_points[0, :] * l0 + y_points[1, :] * l1 + y_points[2, :] * l2
#     return np.dot(y_points.T, all_ls)


