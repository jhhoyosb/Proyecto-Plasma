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


def lagrange(x_points, y_points, x_eval):
    # Generalizar despues...
    l0 = ((x_eval - x_points[1]) / (x_points[0] - x_points[1])) * \
        ((x_eval - x_points[2]) / (x_points[0] - x_points[2]))
    l1 = ((x_eval - x_points[0]) / (x_points[1] - x_points[0])) * \
        ((x_eval - x_points[2]) / (x_points[1] - x_points[2]))
    l2 = ((x_eval - x_points[0]) / (x_points[2] - x_points[0])) * \
        ((x_eval - x_points[1]) / (x_points[2] - x_points[1]))
    return y_points[0, :] * l0 + y_points[1, :] * l1 + y_points[2, :] * l2
