"""Student assignment implementation file.

Complete the TODOs in this file.
"""

# --- Imports --- #
import numpy as np
import scipy.optimize as spo


# --- Problem 1: Quadratic (Projectile Fall Time) --- #
def root_quadratic_projectile_time(h0: float, v0: float, g: float) -> float:
    """Find the positive time at which a vertically launched projectile returns to height 0.

    Parameters
    ----------
    h0 : float
        Initial height.
    v0 : float
        Initial velocity.
    g : float
        Gravitational acceleration.

    Returns
    -------
    float
        Time at which projectile returns to height 0.

    Raises
    ------
    ValueError
        If 'h0' is not greater than or equal to zero.
    ValueError
        If 'g' is zero.
    """

    if h0 < 0:
        raise ValueError("Initial height must be at or above ground level (h=0).")
    elif g == 0:
        raise ValueError("Without gravity, your ball will never touch the ground again...")
    else:
        g = abs(g)
        g = -g * 0.5
        land_t = (-v0 - np.sqrt(v0**2 - 4 * g * h0)) / (2 * g)
        return float(land_t)


# --- Problem 2: Transcendental (Wien's Displacement Law) --- #
def root_wien_displacement_x(coef: float) -> float:
    """Find the positive root x of x * exp(x) = coef * (exp(x) - 1).

    Parameters
    ----------
    coef : float
        Coefficient from the blackbody radiation equation.

    Returns
    -------
    np.float
        Wein's displacement law constant.
    """

    sol = spo.root_scalar(lambda x: x * np.exp(x) - coef * (np.exp(x) - 1), x0=coef)
    return sol.root


# --- Problem 3: Transcendental (Radioactive Decay) --- #
def root_radioactive_decay_time(N0: float, N_target: float, half_life: float) -> float:
    """Find the elapsed time at which a sample of N0 nuclei decays to N_target nuclei.

    Parameters
    ----------
    N0 : float
        Initial sample size.
    N_target : float
        Target sample size.
    half_life : float
        Half-life of sample.

    Returns
    -------
    float
        Elapsed time.
    """

    sol = spo.fsolve(lambda x: N0*(1/2)**(x/half_life) - N_target, half_life)
    return float(sol[0])


# --- Problem 4: Non-Trivial Transcendental (Finite Square Well) --- #
def root_finite_square_well_even_state(V0: float, m: float, hbar: float, L: float) -> float:
    """Find the bound-state wavenumber k of the lowest even state in a finite square well.

    Parameters
    ----------
    V0 : float
        Well depth.
    m : float
        Mass of particle.
    hbar : float
        Value of Dirac constant.
    L : float
        Well width.

    Returns
    -------
    float
        Lowest even bound state's wavenumber.
    """

    def fun(k):
        return np.sqrt((2 * m * V0 / hbar**2) - k**2) - (k * np.tan(k * L * 0.5))
    sol = spo.root(fun, 3 / L)
    return float(sol.x[0])


# --- Problem 5: Complex Roots (Damped Harmonic Oscillator) --- #
def roots_damped_oscillator(m: float, c: float, k: float) -> tuple[complex, complex]:
    """Find both roots s of the characteristic equation m*s**2 + c*s + k = 0.

    Parameters
    ----------
    m : float
        Mass.
    c : float
        Damping factor.
    k : float
        Spring constant.

    Returns
    -------
    tuple(complex, complex)
        Roots of characteristic equation.
    """

    # determine nature of roots using quadratic discriminant
    qdis = c**2 - 4 * m * k
    if qdis > 0:
        root_1 = complex((-c - np.sqrt(qdis)) / (2 * m))
        root_2 = complex((-c + np.sqrt(qdis)) / (2 * m))
        return root_1, root_2
    elif qdis == 0:
        root_d = complex(-c / (2*m))
        return root_d, root_d
    else:
        root_1 = complex((-c - np.emath.sqrt(qdis)) / (2 * m))
        root_2 = complex((-c + np.emath.sqrt(qdis)) / (2 * m))
        return root_1, root_2


# --- Problem 6: Multivariate (Circle Intersection / Trilateration) --- #
def root_circle_intersection(r1: float, r2: float, d: float) -> tuple[float, float]:
    """Find the (x, y) point above the baseline where two circles of radii r1, r2 intersect.

    Parameters
    ----------
    r1 : float
        Radius of first circle.
    r2 : float
        Radius of second circle.
    d : float
        X-offset of second circle.

    Returns
    -------
    tuple[float, float]
        Point above baseline where the two circles meet.
    """

    def fun(x):
        return np.sqrt(r1**2 - x**2) - np.sqrt(r2**2 - (x-d)**2)
    br = [d-r2, r1]
    sol = spo.root_scalar(fun, bracket=br)
    y_val = np.sqrt(r1**2 - (sol.root)**2)
    return tuple([float(sol.root), float(y_val)])


# --- Problem 7: Multivariate Statics (Two Masses, Three Strings) --- #
def root_two_masses_three_strings(
    L: float, lengths: tuple[float, float, float], weights: tuple[float, float]
) -> tuple[tuple[float, float, float], tuple[float, float, float]]:
    """Find the equilibrium angles and tensions of two weights hung from three strings.

    Parameters
    ----------
    L : float
        Horizontal distance between supports.
    lengths : tuple[float, float, float]
        String lengths.
    weights : tuple[float, float]
        Mass of the weights.

    Returns
    -------
    tuple[tuple[float, float, float], tuple[float, float, float]]
        Calculated equilibrium angles and tensions.
    """

    raise NotImplementedError
