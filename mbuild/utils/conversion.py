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


def RB_to_CHARMM(c0, c1, c2, c3, c4, c5):
    """Converts Ryckaert-Bellemans (RB) type dihedrals to CHARMM type
    or

    RB_torsions = c0 + c1*Cos[Psi] + c2*Cos[Psi]^2 + c3*Cos[Psi]^3 + c4*Cos[Psi]^4 + c5*Cos[5*Psi]^5

    where Psi= t-Pi = t - 180 degrees

    Parameters
    ----------
    c0, c1, c2, c3, c4, c5 : Ryckaert-Belleman coefficients (in kcal/mol)

    where:
    n0 = 0
    n1 = 1
    n2 = 2
    n3 = 3
    n4 = 4
    n5 = 5

    and

    d0 = 90
    d1 = 180
    d2 =  0
    d3 =  180
    d4 = 0
    d5 = 180

    converts to:

    CHARMM_torsions =
    = K0 * (1 + Cos[n0*(t) - (d0)] ) + K1 * (1 + Cos[n1*(t) - (d1)] ) + K2 * (1 + Cos[n2*(t) - (d2)] )
    + K3 * (1 + Cos[n3*(t) - (d3)] )  +  K4 * (1 + Cos[n4*(t) - (d4)] )  + K5 * (1 + Cos[n5*(t) - (d5)] )  .

    = K0 + K1 * (1 + Cos[n1*(t) - (d1)] ) + K2 * (1 + Cos[n2*(t) - (d2)] )
    + K3 * (1 + Cos[n3*(t) - (d3)] )  +  K4 * (1 + Cos[n4*(t) - (d4)] )  + K5 * (1 + Cos[n5*(t) - (d5)] )  .

    Returns
    -------
    K0, K1, K2, K3, K4, K5, n0, n1, n2, n3, n4, n5, d0, d1, d2, d3, d4, and d5  : Charmm coefficients (in kcal/mol)

    CHARMM_dihedral coeffs : np.matrix, shape=(6,3)
        Array containing the CHARMM dihedral coeffs  [[K0, n0, d0], [K1, n1, d1], [K2, n2, d2], [K3, n3, d3],
        [K4, n4, d4], [K5, n5, d5]]  (in kcal/mol)

    """
    # see below or the long version is,  K0 = (c0 + c2 / 2 + 3 / 8 * c4) - K1 - K2 - K3 - K4 - K5
    K0 = (c0 - c1 - c3 - (c4/4) - c5)
    K1 = (c1 + (3/4) * c3 + (5/8) * c5)
    K2 =  ((1/2) * c2 + (1/2) * c4)
    K3 =  ((1/4) * c3 + (5/16) * c5)
    K4 = ((1/8) * c4)
    K5 = ((1/16) * c5)


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


    return np.array([[K0, n0, d0], [K1, n1, d1], [K2, n2, d2], [K3, n3, d3], [K4, n4, d4], [K5, n5, d5]])





#***********************************************
# Converting base-10 to base-62 functions (start)
#***********************************************

def base10_to_base62_alph_num(base10_No):
    """Converts base 10 to base 62 so pdb/psf files can add may more than
    9999 atoms and 999 residues."""

    '''base10_No = the base-10 number that you want to convert to base-62)'''

    base62_No = 62
    base10_No = int(base10_No)

    whole_no = 1
    remainder = changeDigit_base10_to_base62_alph_num(int(base10_No % base62_No))
    base62_Values =  str(remainder)
    power = 1

    while whole_no != 0:
        whole_no = int( base10_No / base62_No**power)

        if whole_no == base62_No :
            base62_Values = str(0)+base62_Values

        elif (whole_no != 0) and (whole_no > base62_No) :
            base62_Values = str(changeDigit_base10_to_base62_alph_num(int(whole_no % base62_No))) + base62_Values

        elif (whole_no != 0) and (whole_no < base62_No):
            base62_Values = str(changeDigit_base10_to_base62_alph_num(int(whole_no))) + base62_Values


        power =power+1

    return base62_Values

def changeDigit_base10_to_base62_alph_num(current_digit):
    """The supplemental digits for the base10_to_base62_alph_num function,
    which Converts the base 10 to base 62 """

    """current_digit = the current digit for this base.
    (i.e. in base10 it would be the one, ten, hundreds, or thousands places .....)"""
    base62_values = {j: chr(j+55) if j < 36 else chr(j+61) for j in range(10, 62)}
    return base62_values.get(current_digit, current_digit)




#***********************************************
# Converting base-10 to base-62 functions (end)
#***********************************************





#***********************************************
# Converting base-10 to base-16 functions (start)
#***********************************************

def base10_to_base16_alph_num(base10_No):
    """Converts base 10 to base 16 so pdb/psf files can add may more than
    9999 atoms and 999 residues."""

    """base10_No = the base-10 number that you want to convert to base-16)"""
    return hex(int(base10_No))[2:]


#***********************************************
# Converting base-10 to base-16 functions (end)
#***********************************************


#***********************************************
# Converting base-10 to base-52 functions (start)
#***********************************************

def base10_to_base52_alph(base10_No):
    """Converts base 10 to base 52 so pdb/psf files can add may more than
    atom types in the 3 or 4 character limited pdb and psf files"""

    """base10_No = the base-10 number that you want to convert to base-52)"""

    base52_No = 52
    base10_No = int(base10_No)

    whole_no = 1
    remainder = changeDigit_base10_to_base52_alph_num(int(base10_No % base52_No))
    base52_Values =  str(remainder)
    power = 1

    while whole_no != 0:
        whole_no = int(base10_No / base52_No**power)

        if whole_no == base52_No :
            base52_Values = str('A')+base52_Values

        elif (whole_no != 0) and (whole_no > base52_No) :
            base52_Values =  str(changeDigit_base10_to_base52_alph_num(int(whole_no % base52_No)))+ base52_Values
        elif (whole_no != 0) and (whole_no < base52_No):
            base52_Values =str(changeDigit_base10_to_base52_alph_num(int(whole_no)))+ base52_Values


        power =power+1

    return base52_Values

def changeDigit_base10_to_base52_alph_num(current_digit):
    """The supplemental digits for the base10_to_base52_alph_num function,
    which Converts the base 10 to base 52 """

    """current_digit = the current digit for this base.
    (i.e. in base10 it would be the one, ten, hundreds, or thousands places .....)"""
    base52_values = {j: chr(j+65) if j < 26 else chr(j+71) for j in range(0, 52)}

    current_digit = base52_values[current_digit]
    return current_digit
#***********************************************
# Converting base-10 to base-52 functions (end)
#***********************************************


#***********************************************
# Converting base-10 to base-26 functions (start)
#***********************************************
def base10_to_base26_alph(base10_No):
    """Converts base 10 to base 26 so pdb/psf files can add may more than
    9999 atoms and 999 residues."""

    """base10_No = the base-10 number that you want to convert to base-26)"""

    base26_No = 26
    base10_No = int(base10_No)

    whole_no = 1
    remainder = changeDigit_base10_to_base26_alph_num(int(base10_No % base26_No))
    base26_Values =  str(remainder)
    power = 1

    while whole_no != 0:
        whole_no = int(base10_No / base26_No**power)

        if whole_no == base26_No :
            print('dfe')
            base26_Values = str('A')+base26_Values

        elif (whole_no != 0) and (whole_no > base26_No) :
            base26_Values =  str(changeDigit_base10_to_base26_alph_num(int(whole_no % base26_No)))+ base26_Values
        elif (whole_no != 0) and (whole_no < base26_No):
            base26_Values =str(changeDigit_base10_to_base26_alph_num(int(whole_no)))+ base26_Values


        power =power+1

    return base26_Values

def changeDigit_base10_to_base26_alph_num(current_digit):
    """The supplemental digits for the base10_to_base26_alph_num function,
    which Converts the base 10 to base 26 """

    """current_digit = the current digit for this base.
    (i.e. in base10 it would be the one, ten, hundreds, or thousands places .....)"""
    base26_values = {j: chr(j+65) for j in range(0, 26)}
    current_digit = base26_values[current_digit]
    return current_digit

#***********************************************
# Converting base-10 to base-26 functions (end)
#***********************************************
