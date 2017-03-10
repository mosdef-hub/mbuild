import numpy as np


def RB_to_OPLS(c0, c1, c2, c3, c4, c5):
    """Converts Ryckaert-Bellemans type dihedrals to OPLS type.

    Parameters
    ----------
    c0, c1, c2, c3, c4, c5 : Ryckaert-Belleman coefficients (in kcal/mol)

    Returns
    -------
    opls_coeffs : np.array, shape=(4,)
        Array containing the OPLS dihedrals coeffs f1, f2, f3, and f4
        (in kcal/mol)

    """

    f1 = (-1.5 * c3) - (2 * c1)
    f2 = c0 + c1 + c3
    f3 = -0.5 * c3
    f4 = -0.25 * c4
    return np.array([f1, f2, f3, f4])
