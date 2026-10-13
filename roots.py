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
def root_radioactive_decay_time(N0, N_target, half_life):
    """Find the elapsed time at which a sample of N0 nuclei decays to N_target nuclei."""
    raise NotImplementedError


# --- Problem 4: Non-Trivial Transcendental (Finite Square Well) --- #
def root_finite_square_well_even_state(V0, m, hbar, L):
    """Find the bound-state wavenumber k of the lowest even state in a finite square well."""
    raise NotImplementedError


# --- Problem 5: Complex Roots (Damped Harmonic Oscillator) --- #
def roots_damped_oscillator(m, c, k):
    """Find both roots s of the characteristic equation m*s**2 + c*s + k = 0."""
    raise NotImplementedError


# --- Problem 6: Multivariate (Circle Intersection / Trilateration) --- #
def root_circle_intersection(r1, r2, d):
    """Find the (x, y) point above the baseline where two circles of radii r1, r2 intersect."""
    raise NotImplementedError


# --- Problem 7: Multivariate Statics (Two Masses, Three Strings) --- #
def root_two_masses_three_strings(
    L, lengths, weights
):
    """Find the equilibrium angles and tensions of two weights hung from three strings."""
    raise NotImplementedError
