"""mBuild conversion utilities."""
from warnings import warn

import numpy as np


def RB_to_OPLS(c0, c1, c2, c3, c4, c5, test_mode=False):
    r"""Convert Ryckaert-Bellemans type dihedrals to OPLS type.

    .. math::
    RB_torsions &= c0 + c1*cos(psi) + c2*cos(psi)^2 + c3*cos(psi)^3 + \\
                &= c4*cos(psi)^4 + c5*cos(5*psi)^5

    .. math::
    OPLS_torsions &= f1/2*(1+cos(t)) + f2/2*(1-cos(2*t)) + \\
                  &= f3/2*(1+cos(3*t)) + f4/2(1-cos(4*t))

    where :math:`psi = t - pi = t - 180 degrees`

    Parameters
    ----------
    c0, c1, c2, c3, c4, c5 : Ryckaert-Belleman coefficients (in kcal/mol)
    test_mode : bool, default=False
        Test mode is used to evaluate existing dihedral xml files without
        ending the loops with a 'ValueError' if this function will not work
        for a specific dihedral. This option is needed since every dihedral
        is not guaranteed to convert to OPLS due to the OPLS dihedral lacking
        a 5th power cos term in the OPLS dihedral and lacking a function
        constant (f0 = 0 or f0 must always be zero).
        test_mode SHOULD BE FALSE FOR ALL CASES, EXCEPT TESTING THE XML  FILE FOR ACCURACY.
        False = The function will fail with a ValueError if the RB-torsions
        5th power cos term is not zero (c5 != 0) and the OPLS constant term
        is not zero (f0 != 0).
        True = The function will fail with a ValueError if the RB-torsions
        5th power cos term is not zero (c5 != 0) and the OPLS constant term
        is not zero (f0 != 0).

    Returns
    -------
    opls_coeffs : np.array, shape=(5,)
        Array containing the OPLS dihedrals coeffs f0, f1, f2, f3, and f4
        (in kcal/mol).
        NOTE: fO IS NOT IN THE OPLS DIHEDRAL EQUATION AND IS ONLY USED
        TO TEST IF THIS FUNCTION CAN BE UTILIZED. ONE OF THE
        CONDITIONS OF CONVERSION IS f0 MUST BE ZERO (f0=0).

    Notes
    -----
    c5 must equal zero, or this conversion is not possible.

    (c0 + c1 + c2 + c3 + c4 + c5) must equal zero, or this conversion is not possible.
    """
    if test_mode not in [True, False]:
        raise ValueError(
            "ERROR: test_mode is {}, which is not bool. "
            "Please enter a bool (True or False) for the RB_to_OPLS "
            "functions test_mode option."
            "".format(str(test_mode))
        )

    f0 = 2.0 * (c0 + c1 + c2 + c3 + c4 + c5)

    if bool(np.isclose(c5, 0, atol=1e-12, rtol=0)) is False:
        text_c5_not_zero = (
            "c5 must equal zero, so this conversion is not possible."
        )
        if test_mode is False:
            raise ValueError("ERROR: {}".format(text_c5_not_zero))
        if test_mode is True:
            warn(
                "WARNING: {}"
                "However, the conversion will not fail, and is a warning only "
                "because it is in tesing mode.".format(text_c5_not_zero)
            )

    if bool(np.isclose(f0, 0, atol=1e-12, rtol=0)) is False:
        text_f0_not_zero = (
            "f0 = 2 * (c0 + c1 + c2 + c3 + c4 + c5) must equal zero, "
            "so this conversion is not possible."
        )
        if test_mode is False:
            raise ValueError("ERROR: {}".format(text_f0_not_zero))
        if test_mode is True:
            warn(
                "WARNING: {}"
                "However, the conversion will not fail, and is a warning only "
                "because it is in tesing mode.".format(text_f0_not_zero)
            )

    f1 = -2 * c1 - (3 * c3) / 2
    f2 = -c2 - c4
    f3 = -c3 / 2
    f4 = -c4 / 4
    return np.array([f0, f1, f2, f3, f4])


def OPLS_to_RB(f1, f2, f3, f4):
    r"""Convert OPLS type to Ryckaert-Bellemans type dihedrals.

    .. math::
    OPLS_torsions &= f1/2*(1+cos(t)) + f2/2*(1-cos(2*t)) + \\
                  &= f3/2*(1+cos(3*t)) + f4/2(1-cos(4*t))

    .. math::
    RB_torsions &= c0 + c1*cos(psi) + c2*cos(psi)^2 + c3*cos(psi)^3 + \\
                &= c4*cos(psi)^4 + c5*cos(5*psi)^5

    where :math:`psi = t - pi = t - 180 degrees`

    Parameters
    ----------
    f1, f2, f3, f4 : OPLS dihedrals coeffs (in kcal/mol)

    Returns
    -------
    RB_coeffs : np.array, shape=(6,)
        Array containing the Ryckaert-Bellemans dihedrals
        coeffs c0, c1, c2, c3, c4, and c5 (in kcal/mol)

    """
    c0 = f2 + (f1 + f3) / 2
    c1 = (-f1 + 3 * f3) / 2
    c2 = -f2 + 4 * f4
    c3 = -2 * f3
    c4 = -4 * f4
    c5 = 0
    return np.array([c0, c1, c2, c3, c4, c5])


def RB_to_CHARMM(c0, c1, c2, c3, c4, c5):
    r"""Convert Ryckaert-Bellemans (RB) type dihedrals to CHARMM type.

    .. math::
        RB_torsions &= c0 + c1*cos(psi) + c2*cos(psi)^2 + c3*cos(psi)^3 + \\
                    &= c4*cos(psi)^4 + c5*cos(5*psi)^5

    where :math:`psi = t - pi = t - 180 degrees`

    .. math::
        CHARMM_torsions &= K0 * (1 + cos(n0*t - d0)) + \\
                        &= K1 * (1 + cos(n1*t - d1)) + \\
                        &= K2 * (1 + cos(n2*t - d2)) + \\
                        &= K3 * (1 + cos(n3*t - d3)) + \\
                        &= K4 * (1 + cos(n4*t - d4)) + \\
                        &= K5 * (1 + cos(n5*t - d5))

        CHARMM_torsions &= K0 +
                        &= K1 * (1 + cos(n1*t - d1)) + \\
                        &= K2 * (1 + cos(n2*t - d2)) + \\
                        &= K3 * (1 + cos(n3*t - d3)) + \\
                        &= K4 * (1 + cos(n4*t - d4)) + \\
                        &= K5 * (1 + cos(n5*t - d5))

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
        num = int(number / base ** power)

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
