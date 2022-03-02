"""mBuild conversion utilities."""
from warnings import warn

import numpy as np


def RB_to_OPLS(
    c0,
    c1,
    c2,
    c3,
    c4,
    c5,
    error_tolerance=1e-4,
    error_if_outside_tolerance=True,
):
    r"""Convert Ryckaert-Bellemans type dihedrals to OPLS type.

    .. math::
    RB_{torsions} &= c_0 + c_1*cos(psi) + c_2*cos(psi)^2 + c_3*cos(psi)^3 + \\
                  &= c_4*cos(psi)^4 + c_5*cos(psi)^5

    .. math::
    OPLS_torsions &= \frac{f_0}{2} + \frac{f_1}{2}*(1+cos(t)) + \frac{f_2}{2}*(1-cos(2*t)) + \\
                  &= \frac{f_3}{2}*(1+cos(3*t)) + \frac{f_4}{2}*(1-cos(4*t))

    where :math:`psi = t - pi = t - 180 degrees`

    Parameters
    ----------
    c0, c1, c2, c3, c4, c5 : Ryckaert-Belleman coefficients (in kcal/mol)
    error_tolerance : float, default=1e-4
        The acceptable absolute tolerance between the RB to OPLS conversion.
        Any value entered is converted to an absolute value.
    error_if_outside_tolerance : bool, default=True
        This variable determines whether to provide a ValueError if the RB to OPLS
        conversion is greater than the error_tolerance term (i.e., error_if_outside_tolerance=True),
        or a warning if the RB to OPLS conversion is greater than the error_tolerance term
        (i.e., error_if_outside_tolerance=False).

    Returns
    -------
    opls_coeffs : np.array, shape=(5,)
        Array containing the OPLS dihedrals coeffs f0, f1, f2, f3, and f4
        (in kcal/mol).


    Notes
    -----
    c5 must equal zero, or this conversion is not possible.

    (c0 + c1 + c2 + c3 + c4 + c5) must equal zero, to have the exact
    energies represented in the dihedral.

    NOTE: fO IS TYPICALLY NOT IN THE OPLS DIHEDRAL EQUATION AND IS
    ONLY USED TO TEST IF THIS FUNCTION CAN BE UTILIZED, AND
    DETERMINE THE fO VALUE FOR LATER ENERGY SCALING IN MOLECULAR DYNAMICS
    (MD) SIMULATIONS. THIS FUNCTION TESTS IF f0 IS ZERO (f0=0).

    .. warning:: The :math:`\frac{f_{0}}{2}` term is the constant for the OPLS dihedral equation.
        If the f0 term is not zero, the dihedral is not an exact conversion;
        since this constant does not contribute to the force equation,
        this should provide matching results for MD, but the energy for each
        dihedral will be shifted by the :math:`\frac{f_{0}}{2}` value.
    """
    if not isinstance(error_tolerance, float):
        raise TypeError(
            f"The error_tolerance variable must be a float, is type {type(error_tolerance)}."
        )
    error_tolerance = abs(error_tolerance)

    if not isinstance(error_if_outside_tolerance, bool):
        raise TypeError(
            f"The text_for_error_tolerance variable must be a bool, is type {type(error_if_outside_tolerance)}."
        )

    if not np.all(np.isclose(c5, 0, atol=1e-10, rtol=0)):
        raise ValueError(
            "c5 must equal zero, so this conversion is not possible."
        )

    f0 = 2.0 * (c0 + c1 + c2 + c3 + c4 + c5)
    if not np.all(np.isclose(f0 / 2, 0, atol=error_tolerance, rtol=0)):
        text_for_error_tolerance = (
            "f0 = 2 * ( c0 + c1 + c2 + c3 + c4 + c5 ) is not zero. "
            "The f0/2 term is the constant for the OPLS dihedral. "
            "Since the f0 term is not zero, the dihedral is not an "
            "exact conversion; since this constant does not contribute "
            "to the force equation, this should provide matching results "
            "for MD, but the energy for each dihedral will be shifted "
            "by the f0/2 value."
        )
        if error_if_outside_tolerance is True:
            raise ValueError(text_for_error_tolerance)
        elif error_if_outside_tolerance is False:
            warn(text_for_error_tolerance)

    f1 = -2 * c1 - (3 * c3) / 2
    f2 = -c2 - c4
    f3 = -c3 / 2
    f4 = -c4 / 4
    return np.array([f0, f1, f2, f3, f4])


def OPLS_to_RB(f0, f1, f2, f3, f4, error_tolerance=1e-4):
    r"""Convert OPLS type to Ryckaert-Bellemans type dihedrals.

    .. math::
    OPLS_torsions &= \frac{f_0}{2} + \frac{f_1}{2}*(1+cos(t)) + \frac{f_2}{2}*(1-cos(2*t)) + \\
                  &= \frac{f_3}{2}*(1+cos(3*t)) + \frac{f_4}{2}*(1-cos(4*t))

    .. math::
    RB_{torsions} &= c_0 + c_1*cos(psi) + c_2*cos(psi)^2 + c_3*cos(psi)^3 + \\
                  &= c_4*cos(psi)^4 + c_5*cos(psi)^5

    where :math:`psi = t - pi = t - 180 degrees`

    Parameters
    ----------
    f0, f1, f2, f3, f4 : OPLS dihedrals coeffs (in kcal/mol)
    error_tolerance : float, default=1e-4
        The acceptable absolute tolerance between the OPLS to RB conversion
        without throwing a warning. Any value entered is converted to an
        absolute value.

    Returns
    -------
    RB_coeffs : np.array, shape=(6,)
        Array containing the Ryckaert-Bellemans dihedrals
        coeffs c0, c1, c2, c3, c4, and c5 (in kcal/mol)

    Notes
    -----
    NOTE: fO IS TYPICALLY NOT IN THE OPLS DIHEDRAL EQUATION (i.e., f0=0).

    .. warning:: The :math:`\frac{f_{0}}{2}` term is the constant for the OPLS dihedral equation,
        which is and added to a constant for the RB torsions equation via the c0 coefficient.
        If the f0 term is zero in the OPLS dihedral form or is force set to zero in this equation,
        the dihedral is may not an exact conversion;
        since this constant does not contribute to the force equation,
        this should provide matching results for MD, but the energy for each
        dihedral will be shifted by the real :math:`\frac{f_{0}}{2}` value.
    """
    if not isinstance(error_tolerance, float):
        raise TypeError(
            f"The error_tolerance variable must be a float, is type {type(error_tolerance)}."
        )
    error_tolerance = abs(error_tolerance)

    if np.all(np.isclose(f0 / 2, 0, atol=error_tolerance, rtol=0)):
        warn(
            "The f0/2 term is the constant for the OPLS dihedral equation, "
            "which is added to a constant for the RB torsions equation via the c0 coefficient. "
            "The f0 term is zero in the OPLS dihedral form or is force set to zero in this equation, "
            "so the dihedral is may not an exact conversion; "
            "since this constant does not contribute to the force equation, "
            "this should provide matching results for MD, but the energy for each"
            "dihedral will be shifted by the real f0/2 value."
        )

    c0 = f2 + (f0 + f1 + f3) / 2
    c1 = (-f1 + 3 * f3) / 2
    c2 = -f2 + 4 * f4
    c3 = -2 * f3
    c4 = -4 * f4
    c5 = 0
    return np.array([c0, c1, c2, c3, c4, c5])


def RB_to_CHARMM(c0, c1, c2, c3, c4, c5):
    r"""Convert Ryckaert-Bellemans (RB) type dihedrals to CHARMM type.

    .. math::
        RB_torsions &= c_0 + c_1*cos(psi) + c_2*cos(psi)^2 + c_3*cos(psi)^3 + \\
                    &= c_4*cos(psi)^4 + c_5*cos(psi)^5

    where :math:`psi = t - pi = t - 180 degrees`

    .. math::
        CHARMM_torsions &= K_0 * (1 + cos(n_0*t - d_0)) + \\
                        &= K_1 * (1 + cos(n_1*t - d_1)) + \\
                        &= K_2 * (1 + cos(n_2*t - d_2)) + \\
                        &= K_3 * (1 + cos(n_3*t - d_3)) + \\
                        &= K_4 * (1 + cos(n_4*t - d_4)) + \\
                        &= K_5 * (1 + cos(n_5*t - d_5))

        CHARMM_torsions &= K_0 +
                        &= K_1 * (1 + cos(n_1*t - d_1)) + \\
                        &= K_2 * (1 + cos(n_2*t - d_2)) + \\
                        &= K_3 * (1 + cos(n_3*t - d_3)) + \\
                        &= K_4 * (1 + cos(n_4*t - d_4)) + \\
                        &= K_5 * (1 + cos(n_5*t - d_5))

    Parameters
    ----------
    c0, c1, c2, c3, c4, c5 : Ryckaert-Belleman coefficients (in kcal/mol)

    n0 = 0
    n1 = 1
    n2 = 2
    n3 = 3
    n4 = 4
    n5 = 5

    d0 = 90
    d1 = 180
    d2 = 0
    d3 = 180
    d4 = 0
    d5 = 180

    Returns
    -------
    CHARMM_dihedral coeffs : np.matrix, shape=(6,3)
        Array containing the CHARMM dihedral coeffs (in kcal/mol):

        [[K0, n0, d0],
         [K1, n1, d1],
         [K2, n2, d2],
         [K3, n3, d3],
         [K4, n4, d4],
         [K5, n5, d5]]
    """
    # see below or the long version is,
    # K0 = (c0 + c2 / 2 + 3 / 8 * c4) - K1 - K2 - K3 - K4 - K5
    K0 = c0 - c1 - c3 - (c4 / 4) - c5
    K1 = c1 + (3 / 4) * c3 + (5 / 8) * c5
    K2 = (1 / 2) * c2 + (1 / 2) * c4
    K3 = (1 / 4) * c3 + (5 / 16) * c5
    K4 = (1 / 8) * c4
    K5 = (1 / 16) * c5

    n0 = 0
    n1 = 1
    n2 = 2
    n3 = 3
    n4 = 4
    n5 = 5

    d0 = 90
    d1 = 180
    d2 = 0
    d3 = 180
    d4 = 0
    d5 = 180

    return np.array(
        [
            [K0, n0, d0],
            [K1, n1, d1],
            [K2, n2, d2],
            [K3, n3, d3],
            [K4, n4, d4],
            [K5, n5, d5],
        ]
    )


def base10_to_base62_alph_num(base10_no):
    """Convert base-10 integer to base-62 alphanumeric system.

    This function provides a utility to write pdb/psf files such that it can
    add may more than 9999 atoms and 999 residues.

    Parameters
    ----------
    base10_no: int
        The integer to convert to base-62 alphanumeric system

    Returns
    -------
    str
        The converted base-62 system string

    See Also
    --------
    mbuild.conversion._to_base: Helper function to perform a base-n conversion
    """
    return _to_base(base10_no, base=62)


def base10_to_base52_alph(base10_no):
    """Convert base-10 integer to base-52 alphabetic system.

    This function provides a utility to write pdb/psf files such that it can
    add more atom types in the 3 or 4 character limited pdb and psf files

    Parameters
    ----------
    base10_no: int
        The integer to convert to base-52 alphabetic system

    Returns
    -------
    str
        The converted base-52 system string

    See Also
    --------
    mbuild.conversion._to_base: Helper function to perform a base-n conversion
    """
    return _to_base(number=base10_no, base=52)


def base10_to_base26_alph(base10_no):
    """Convert base-10 integer to base-26 alphabetic system.

    This function provides a utility to write pdb/psf files such that it can
    add many more than 9999 atoms and 999 residues.

    Parameters
    ----------
    base10_no: int
        The integer to convert to base-26 alphabetic system

    Returns
    -------
    str
        The converted base-26 system string

    See Also
    --------
    mbuild.conversion._to_base: Helper function to perform a base-n conversion
    """
    return _to_base(base10_no, base=26)


def base10_to_base16_alph_num(base10_no):
    """Convert base-10 integer to base-16 hexadecimal system.

    This function provides a utility to write pdb/psf files such that it can
    add many more than 9999 atoms and 999 residues.

    Parameters
    ----------
    base10_no: int
        The integer to convert to base-16 hexadecimal system

    Returns
    -------
    str
        The converted base-16 system string

    See Also
    --------
    mbuild.conversion._to_base: Helper function to perform a base-n conversion
    """
    return hex(int(base10_no))[2:]


# Helpers to convert base
def _to_base(number, base=62):
    """Convert a base-10 number into base-n alpha-num."""
    start_values = {62: "0", 52: "A", 26: "A"}
    if base not in start_values:
        raise ValueError(
            f"Base-{base} system is not supported. Supported bases are: "
            f"{list(start_values.keys())}"
        )

    num = 1
    number = int(number)
    remainder = _digit_to_alpha_num((number % base), base)
    base_n_values = str(remainder)
    power = 1

    while num != 0:
        num = int(number / base**power)

        if num == base:
            base_n_values = start_values[base] + base_n_values

        elif num != 0 and num > base:
            base_n_values = (
                str(_digit_to_alpha_num(int(num % base), base)) + base_n_values
            )

        elif (num != 0) and (num < base):
            base_n_values = (
                str(_digit_to_alpha_num(int(num), base)) + base_n_values
            )

        power += 1

    return base_n_values


def _digit_to_alpha_num(digit, base=52):
    """Convert digit to base-n."""
    base_values = {
        26: {j: chr(j + 65) for j in range(0, 26)},
        52: {j: chr(j + 65) if j < 26 else chr(j + 71) for j in range(0, 52)},
        62: {j: chr(j + 55) if j < 36 else chr(j + 61) for j in range(10, 62)},
    }

    if base not in base_values:
        raise ValueError(
            f"Base-{base} system is not supported. Supported bases are: "
            f"{list(base_values.keys())}"
        )

    return base_values[base].get(digit, digit)
