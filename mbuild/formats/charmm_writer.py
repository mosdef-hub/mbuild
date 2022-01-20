import datetime
import os
from collections import OrderedDict
from warnings import warn

import numpy as np
from parmed.periodic_table import Element
from parmed.utils.io import genopen

from mbuild.box import Box
from mbuild.compound import Compound
from mbuild.utils.conversion import (
    RB_to_CHARMM,
    base10_to_base16_alph_num,
    base10_to_base26_alph,
    base10_to_base52_alph,
    base10_to_base62_alph_num,
)
from mbuild.utils.sorting import natural_sort
from mbuild.utils.specific_ff_to_residue import specific_ff_to_residue


def _get_bond_type_key(
    bond, sigma_conversion_factor, epsilon_conversion_factor
):
    """Get the bond_type key for a bond

    Parameters
    ----------
    bond : mbuild.compound.Compound
        The bond information from the mbuild.compound.Compound
    sigma_conversion_factor : float or int
        The sigma conversion factor
    epsilon_conversion_factor : float or int
        The epsilon conversion factor

    Returns
    ----------
    tuple, (bond_k_constant, bond_bo_length, bond_atom_1_and_2_types_tuple,
        bond_atom_1_residue_name, bond_atom_2_residue_name)
        bond_k_constant : float
            Harmonic bond k-constant or bond energy scaling factor
        bond_bo_length : float
            Harmonic bonds equilibrium length
        bond_atom_1_and_2_types_tuple : tuple
            A sorted tuple ofstrings for the bonded atom types for atoms 1 and 2.
        bond_atom_1_residue_name : str
            The residue name for atom 1 in the bond.
        bond_atom_2_residue_name : str
            The residue name for atom 2 in the bond.
    """
    bond_k_constant = round(
        bond.type.k
        * (sigma_conversion_factor ** 2 / epsilon_conversion_factor),
        8,
    )
    bond_bo_length = round(bond.type.req / sigma_conversion_factor, 8)
    bond_atom_1_and_2_types_tuple = tuple(
        sorted((bond.atom1.type, bond.atom2.type))
    )
    bond_atom_1_residue_name = bond.atom1.residue.name
    bond_atom_2_residue_name = bond.atom2.residue.name

    return (
        bond_k_constant,
        bond_bo_length,
        bond_atom_1_and_2_types_tuple,
        bond_atom_1_residue_name,
        bond_atom_2_residue_name,
    )


def _get_angle_type_key(
    angle, sigma_conversion_factor, epsilon_conversion_factor
):
    """Get the angle_type key for an harmonic angle

    Parameters
    ----------
    angle : parmed.topologyobjects.Angle
        The angle information from the parmed.topologyobjects.Angle
    sigma_conversion_factor : float or int
        The sigma conversion factor
    epsilon_conversion_factor : float or int
        The epsilon conversion factor

    Returns
    ----------
    tuple, (angle_k_constant, angle_theta_o, angle_center_atom_type_2, angle_end_atom_types_1_and_3_tuple,
        angle_residue_atom_1, angle_residue_atom_2, angle_residue_atom_3)
        angle_k_constant : float
            Harmonic angle k-constant or bond energy scaling factor
        angle_theta_o : float
            Harmonic equilbrium angle between the atoms
        angle_center_atom_type_2 : str
            The center atom type for the angle (atom type 2)
        angle_end_atom_types_1_and_3_tuple : tuple
            A sorted tuple of atom types (strings) for the end angle atoms
            (atoms 1 and 3).
        angle_atom_1_residue_name : str
            The residue name for atom 1 in the angle (end atom).
        angle_atom_2_residue_name : str
            The residue name for atom 2 in the angle (center atom).
        angle_atom_3_residue_name : str
            The residue name for atom 3 in the angle (end atom).
    """

    angle_k_constant = round(
        angle.type.k
        * (sigma_conversion_factor ** 2 / epsilon_conversion_factor),
        8,
    )
    angle_theta_o = round(angle.type.theteq, 8)
    angle_center_atom_type_2 = angle.atom2.type
    angle_end_atom_types_1_and_3_tuple = tuple(
        sorted((angle.atom1.type, angle.atom3.type))
    )
    angle_residue_atom_1 = angle.atom1.residue.name
    angle_residue_atom_2 = angle.atom2.residue.name
    angle_residue_atom_3 = angle.atom3.residue.name

    return (
        angle_k_constant,
        angle_theta_o,
        angle_center_atom_type_2,
        angle_end_atom_types_1_and_3_tuple,
        angle_residue_atom_1,
        angle_residue_atom_2,
        angle_residue_atom_3,
    )


def _get_dihedral_rb_torsion_key(dihedral, epsilon_conversion_factor):
    """Get the dihedral_type key for a Ryckaert-Bellemans (RB) dihedrals/torsions

    Parameters
    ----------
    dihedral : parmed.topologyobjects.Dihedral
        The dihedral information from the parmed.topologyobjects.Angle
    epsilon_conversion_factor : float or int
        The epsilon conversion factor

    Returns
    ----------
    tuple, (dihed_type_RB_c0,  dihed_type_RB_c1, dihed_type_RB_c2, dihed_type_RB_c3, dihed_type_RB_c4,
        dihed_type_RB_c5, dihed_type_scee, dihed_type_scnb, dihed_atom_1_type, dihed_atom_2_type,
        dihed_atom_3_type, dihed_atom_4_type, dihed_atom_1_res_type, dihed_atom_2_res_type,
        dihed_atom_3_res_type,  dihed_atom_4_res_type)
        dihed_type_RB_c0 : float
            Ryckaert-Bellemans (RB) dihedrals/torsions C0 constant.
        dihed_type_RB_c1 : float
            Ryckaert-Bellemans (RB) dihedrals/torsions C1 constant.
        dihed_type_RB_c2 : float
            Ryckaert-Bellemans (RB) dihedrals/torsions C2 constant.
        dihed_type_RB_c3 : float
            Ryckaert-Bellemans (RB) dihedrals/torsions C3 constant.
        dihed_type_RB_c4 : float
            Ryckaert-Bellemans (RB) dihedrals/torsions C4 constant.
        dihed_type_RB_c5 : float
            Ryckaert-Bellemans (RB) dihedrals/torsions C5 constant.
        dihed_type_scee : float, default = 1.0
            The 1-4 electrostatic scaling factor
        dihed_type_scnb : float, default = 1.0
            The 1-4 Lennard-Jones scaling factor.
        dihed_atom_1_type : str
            The atom type for atom number 1 in the dihedral
        dihed_atom_2_type : str
            The atom type for atom number 2 in the dihedral
        dihed_atom_3_type : str
            The atom type for atom number 3 in the dihedral
        dihed_atom_4_type : str
            The atom type for atom number 4 in the dihedral
        dihed_atom_1_res_type : str
            The residue name for atom number 1 in the dihedral
        dihed_atom_2_res_type : str
            The residue name for atom number 2 in the dihedral
        dihed_atom_3_res_type : str
            The residue name for atom number 3 in the dihedral
        dihed_atom_4_res_type : str
            The residue name for atom number 4 in the dihedral
    """

    lj_unit = 1 / epsilon_conversion_factor

    dihed_type_RB_c0 = round(dihedral.type.c0 * lj_unit, 8)
    dihed_type_RB_c1 = round(dihedral.type.c1 * lj_unit, 8)
    dihed_type_RB_c2 = round(dihedral.type.c2 * lj_unit, 8)
    dihed_type_RB_c3 = round(dihedral.type.c3 * lj_unit, 8)
    dihed_type_RB_c4 = round(dihedral.type.c4 * lj_unit, 8)
    dihed_type_RB_c5 = round(dihedral.type.c5 * lj_unit, 8)

    dihed_type_scee = round(dihedral.type.scee, 4)
    dihed_type_scnb = round(dihedral.type.scnb, 4)

    dihed_atom_1_type = dihedral.atom1.type
    dihed_atom_2_type = dihedral.atom2.type
    dihed_atom_3_type = dihedral.atom3.type
    dihed_atom_4_type = dihedral.atom4.type

    dihed_atom_1_res_type = dihedral.atom1.residue.name
    dihed_atom_2_res_type = dihedral.atom2.residue.name
    dihed_atom_3_res_type = dihedral.atom3.residue.name
    dihed_atom_4_res_type = dihedral.atom4.residue.name

    return (
        dihed_type_RB_c0,
        dihed_type_RB_c1,
        dihed_type_RB_c2,
        dihed_type_RB_c3,
        dihed_type_RB_c4,
        dihed_type_RB_c5,
        dihed_type_scee,
        dihed_type_scnb,
        dihed_atom_1_type,
        dihed_atom_2_type,
        dihed_atom_3_type,
        dihed_atom_4_type,
        dihed_atom_1_res_type,
        dihed_atom_2_res_type,
        dihed_atom_3_res_type,
        dihed_atom_4_res_type,
    )


def _get_improper_type_key(improper, epsilon_conversion_factor):
    """Get the improper_type key for the harmonic improper

    Parameters
    ----------
    improper : parmed.topologyobjects.Dihedral
        The improper information from the parmed.topologyobjects.Angle
    epsilon_conversion_factor : float or int
        The epsilon conversion factor

    Returns
    ----------
    tuple, (improper_k_constant, improper_psi_o, improper_atom_1_type,
        (improper_atom_2_type, improper_atom_3_type,  improper_atom_4_type), improper_atom_1_res_type,
        (improper_atom_2_res_type, improper_atom_3_res_type, improper_atom_4_res_type)
        improper_k_constant : float
            Harmonic k-constant or bond energy scaling factor
        improper_psi_o : float
            Harmonic equilbrium improper angle
        improper_atom_1_type : str
            The atom type for atom number 1 in the dihedral
        improper_atom_2_type : str
            The atom type for atom number 2 in the dihedral
        improper_atom_3_type : str
            The atom type for atom number 3 in the dihedral
        improper_atom_4_type : str
            The atom type for atom number 4 in the dihedral
        improper_atom_1_res_type : str
            The residue name for atom number 1 in the dihedral
        improper_atom_2_res_type : str
            The residue name for atom number 2 in the dihedral
        improper_atom_3_res_type : str
            The residue name for atom number 3 in the dihedral
        improper_atom_4_res_type : str
            The residue name for atom number 4 in the dihedral
    """
    lj_unit = 1 / epsilon_conversion_factor

    improper_k_constant = round(improper.type.psi_k * lj_unit, 8)
    improper_psi_o = round(improper.type.psi_eq, 8)
    improper_atom_1_type = improper.atom1.type
    improper_atom_2_type = improper.atom2.type
    improper_atom_3_type = improper.atom3.type
    improper_atom_4_type = improper.atom4.type
    improper_atom_1_res_type = improper.atom1.residue.name
    improper_atom_2_res_type = improper.atom2.residue.name
    improper_atom_3_res_type = improper.atom3.residue.name
    improper_atom_4_res_type = improper.atom4.residue.name

    return (
        improper_k_constant,
        improper_psi_o,
        improper_atom_1_type,
        (improper_atom_2_type, improper_atom_3_type, improper_atom_4_type),
        improper_atom_1_res_type,
        (
            improper_atom_2_res_type,
            improper_atom_3_res_type,
            improper_atom_4_res_type,
        ),
    )


def _get_unique_bond_types(
    structure, sigma_conversion_factor, epsilon_conversion_factor
):
    """Get the unique bond types for a structure in a dictionary

    Parameters
    ----------
    structure : parmed.structure.Structure
       This is a parmed stucture input (parmed.structure.Structure)
    sigma_conversion_factor : float or int
        The sigma conversion factor
    epsilon_conversion_factor : float or int
        The epsilon conversion factor

    Returns
    ----------
    bond_key_dict : dict, {(float, float, (str, str), str, str) : unique_number}
        This provides a way to uniquely number the harmonic bond types
        by providing all the bond parameters as a key and the
        unique number as the value. An example of the dict is below:
        {(bond_k_constant, bond_bo_length, (bond_atom_type_1, bond_atom_type_2),
          bond_atom_1_residue_name, bond_atom_2_residue_name ) : 1,
         ...,
         (bond_k_constant, bond_bo_length, (bond_atom_type_1, bond_atom_type_2),
          bond_atom_1_residue_name, bond_atom_2_residue_name ) : n
        }
    """

    unique_bond_set = set()
    for bond in structure.bonds:
        unique_bond_set.add(
            _get_bond_type_key(
                bond, sigma_conversion_factor, epsilon_conversion_factor
            )
        )

    bond_key_dict = {
        bond_key: i + 1 for i, bond_key in enumerate(unique_bond_set)
    }

    return bond_key_dict


def _get_unique_angle_types(
    structure, sigma_conversion_factor, epsilon_conversion_factor
):
    """Get the unique angle types for a structure and return a dictionary

    Parameters
    ----------
    structure : parmed.structure.Structure
       This is a parmed stucture input (parmed.structure.Structure)
    sigma_conversion_factor : float or int
        The sigma conversion factor
    epsilon_conversion_factor : float or int
        The epsilon conversion factor

    Returns
    ----------
    angle_key_dict : dict, {(float, float, str, (str, str), str, str, str) : unique_number}
        This provides a way to uniquely number the harmonic angle types
        by providing all the angle parameters as a key and the
        unique number as the value. An example of the dict is below:
        {(angle_k_constant, angle_theta_o, angle_center_atom_type_2,
          (angle_end_atom_type_1, angle_end_atom_type_3),
           angle_residue_atom_1, angle_residue_atom_2, angle_residue_atom_3
         ) : 1,
         ...,
         (angle_k_constant, angle_theta_o, angle_center_atom_type_2,
         (angle_end_atom_type_1, angle_end_atom_type_3),
          angle_residue_atom_1, angle_residue_atom_2, angle_residue_atom_3) : n
        }
    """

    unique_angle_set = set()
    for angle in structure.angles:
        unique_angle_set.add(
            _get_angle_type_key(
                angle, sigma_conversion_factor, epsilon_conversion_factor
            )
        )

    angle_key_dict = {
        angle_key: i + 1 for i, angle_key in enumerate(unique_angle_set)
    }

    return angle_key_dict


def _get_unique_rb_torsion_types(structure, epsilon_conversion_factor):
    """Get the unique rb torsion types for a structure and return a dictionary

    Parameters
    ----------
    structure : parmed.structure.Structure
       This is a parmed stucture input (parmed.structure.Structure)
    epsilon_conversion_factor : float or int
        The epsilon conversion factor

    Returns
    ----------
    dihed_key_dict : dict, {(float, float, float, float, float, float, float, float,
                             str, str, str, str, str, str, str, str) : unique_number}
        This provides a way to uniquely number the Ryckaert-Bellemans (RB)
        dihedral types by providing all the dihedral parameters as a key and the
        unique number as the value. An example of the dict is below:
        {(dihed_type_RB_c0, dihed_type_RB_c1, dihed_type_RB_c2,
          dihed_type_RB_c3, dihed_type_RB_c4, dihed_type_RB_c5,
          dihed_type_scee, dihed_type_scnb,
          dihed_atom_1_type, dihed_atom_2_type,
          dihed_atom_3_type, dihed_atom_4_type,
          dihed_atom_1_res_type, dihed_atom_2_res_type,
          dihed_atom_3_res_type, dihed_atom_4_res_type
         ) : 1,
         ...,
        (dihed_type_RB_c0, dihed_type_RB_c1, dihed_type_RB_c2,
          dihed_type_RB_c3, dihed_type_RB_c4, dihed_type_RB_c5,
          dihed_type_scee, dihed_type_scnb,
          dihed_atom_1_type, dihed_atom_2_type,
          dihed_atom_3_type, dihed_atom_4_type,
          dihed_atom_1_res_type, dihed_atom_2_res_type,
          dihed_atom_3_res_type, dihed_atom_4_res_type
         ) : n
        }
    """
    unique_dihedral_set = set()
    for dihedral in structure.rb_torsions:
        unique_dihedral_set.add(
            _get_dihedral_rb_torsion_key(dihedral, epsilon_conversion_factor)
        )

    dihed_key_dict = {
        dihed_key: i + 1 for i, dihed_key in enumerate(unique_dihedral_set)
    }

    return dihed_key_dict


def _get_unique_improper_types(structure, epsilon_conversion_factor):
    """Get the unique improper types for a structure  and return a dictionary

    Parameters
    ----------
    structure : parmed.structure.Structure
       This is a parmed stucture input (parmed.structure.Structure)
    epsilon_conversion_factor : float or int
        The epsilon conversion factor

    Returns
    ----------
    improper_key_dict : dict, {(float, float, str, (str, str, str), str, (str, str, str)) : unique_number}
        This provides a way to uniquely number the harmonic improper
        types by providing all the improper parameters as a key and the
        unique number as the value. An example of the dict is below:
        {(improper_k_constant, improper_psi_o,
          improper_atom_1_type, (improper_atom_2_type,
          improper_atom_3_type, improper_atom_4_type),
          improper_atom_1_res_type, (improper_atom_2_res_type,
          improper_atom_3_res_type, improper_atom_4_res_type)
         ) : 1,
         ...,
         (improper_k_constant, improper_psi_o,
          improper_atom_1_type, (improper_atom_2_type,
          improper_atom_3_type, improper_atom_4_type),
          improper_atom_1_res_type, (improper_atom_2_res_type,
          improper_atom_3_res_type, improper_atom_4_res_type)
         ) : n
        }
    """
    unique_improper_set = set()
    for improper in structure.impropers:
        unique_improper_set.add(
            _get_improper_type_key(improper, epsilon_conversion_factor)
        )

    improper_key_dict = {
        improper_key: i + 1
        for i, improper_key in enumerate(unique_improper_set)
    }

    return improper_key_dict


def _get_bond_types(
    structure, sigma_conversion_factor, epsilon_conversion_factor
):
    """Get a list of unique bond types that are used to create the Charmm style
    parameter file (i.e. force field file).  This also checks that the reverse bond order
    is considered the same unique bond type.

    Parameters
    ----------
    structure : parmed.Structure
        The parmed structure
    sigma_conversion_factor : float or int
        The sigma conversion factor
    epsilon_conversion_factor : float or int
        The epsilon conversion factor

    Returns
    ----------
    bond_types : list
        The bond types by number in the structure
    unique_bond_types : OrderedDict, ((float, float, (str, str), str, str), unique_number)
        This provides the unique harmonic bond types, numbering, and the data values
        so it can easily be extracted.  An example of the OrderedDict is below:
        OrderedDict([(bond_k_constant, bond_bo_length, (bond_atom_type_1, bond_atom_type_2),
                      bond_atom_1_residue_name, bond_atom_2_residue_name), 1),
                     ...,
                    (bond_k_constant, bond_bo_length, (bond_atom_type_1, bond_atom_type_2),
                     bond_atom_1_residue_name, bond_atom_2_residue_name), n)])
    """

    unique_bond_types = _get_unique_bond_types(
        structure, sigma_conversion_factor, epsilon_conversion_factor
    )

    bond_types = [
        unique_bond_types[
            _get_bond_type_key(
                bond, sigma_conversion_factor, epsilon_conversion_factor
            )
        ]
        for bond in structure.bonds
    ]

    unique_bond_check_dict = {}
    for i_value_bond, i_key_bond in unique_bond_types.items():
        i_value_duplicated = False
        for j_value_bond, j_key_bond in unique_bond_types.items():
            j_value_bond_reorder = (
                j_value_bond[0],
                j_value_bond[1],
                j_value_bond[2][0],
                j_value_bond[2][0],
                j_value_bond[3],
                j_value_bond[4],
            )

            if i_value_bond == j_value_bond_reorder:
                i_value_duplicated = True
                if i_value_bond[2][0] > j_value_bond[2][0]:
                    unique_bond_check_dict.update(
                        {j_value_bond: len(unique_bond_check_dict)}
                    )
                else:
                    unique_bond_check_dict.update(
                        {i_value_bond: len(unique_bond_check_dict)}
                    )

            if i_value_duplicated is False:
                unique_bond_check_dict.update(
                    {i_value_bond: len(unique_bond_check_dict)}
                )

    unique_bond_types = OrderedDict(
        [(y, x) for y, x in unique_bond_check_dict.items()]
    )

    return bond_types, unique_bond_types


def _get_angle_types(
    structure,
    sigma_conversion_factor,
    epsilon_conversion_factor,
    use_urey_bradleys=False,
):
    """
    Get a list of unique angle types that are used to create the Charmm style
    parameter file (i.e. force field file).  This also checks that the alternately
    ordered angle types are considered the same unique angle type.

    Parameters
    ----------
    structure : parmed.Structure
        The parmed structure
    sigma_conversion_factor : float or int
        The sigma conversion factor
    epsilon_conversion_factor : float or int
        The epsilon conversion factor
    use_urey_bradleys : bool
        The option that Urey-Bradleys are included in the angles

    Returns
    ----------
    angle_types : list
        The angle types by number in the structure
    unique_angle_types : OrderedDict, ((float, float, (str, str), str, str), unique_number)
        This provides the unique harmonic angles types, numbering, and the data values
        so it can easily be extracted.  An example of the OrderedDict is below:
        OrderedDict([(angle_k_constant, angle_theta_o, angle_center_atom_type_2,
                     (angle_end_atom_type_1, angle_end_atom_type_3),
                      angle_residue_atom_1, angle_residue_atom_2, angle_residue_atom_3), 1,
                     ...,
                     (angle_k_constant, angle_theta_o, angle_center_atom_type_2,
                     (angle_end_atom_type_1, angle_end_atom_type_3),
                     angle_residue_atom_1, angle_residue_atom_2, angle_residue_atom_3), n])
    """

    if use_urey_bradleys:
        print_warn_text = (
            "WARNING: Urey-Bradleys are not available in the current "
            "version of this psf, pdb, and GOMC writer."
        )
        warn(print_warn_text)
        return None, None
    else:
        unique_angle_types = _get_unique_angle_types(
            structure, sigma_conversion_factor, epsilon_conversion_factor
        )

        angle_types = [
            unique_angle_types[
                _get_angle_type_key(
                    angle, sigma_conversion_factor, epsilon_conversion_factor
                )
            ]
            for angle in structure.angles
        ]

    unique_angle_check_dict = {}
    for i_value_ang, i_key_ang in unique_angle_types.items():
        i_value_duplicated = False
        for j_value_ang, j_key_ang in unique_angle_types.items():
            j_value_ang_reorder = (
                j_value_ang[0],
                j_value_ang[1],
                j_value_ang[2],
                j_value_ang[3][0],
                j_value_ang[3][1],
                j_value_ang[4],
                j_value_ang[5],
                j_value_ang[6],
            )

            if i_value_ang == j_value_ang_reorder:
                i_value_duplicated = True
                if i_value_ang[2] > j_value_ang[2]:
                    unique_angle_check_dict.update(
                        {j_value_ang: len(unique_angle_check_dict)}
                    )
                else:
                    unique_angle_check_dict.update(
                        {i_value_ang: len(unique_angle_check_dict)}
                    )

        if not i_value_duplicated:
            unique_angle_check_dict.update(
                {i_value_ang: len(unique_angle_check_dict)}
            )

    unique_angle_types = OrderedDict(
        [(y, x) for y, x in unique_angle_check_dict.items()]
    )

    return angle_types, unique_angle_types


def _get_dihedral_types(
    structure, use_rb_torsions, use_dihedrals, epsilon_conversion_factor
):
    """
    Get a list of unique dihedral types that are used to create the Charmm style
    parameter file (i.e. force field file).  This also checks that the alternately
    ordered dihedral types are considered the same unique dihedral type.

    Parameters
    ----------
    structure : parmed.Structure
        The parmed structure
    use_rb_torsions : bool
        The Ryckaert-Bellemans (RB) dihedrals/torsions
    use_dihedrals : bool
        The CHARMM style dihedrals style equations.
    epsilon_conversion_factor : float or int
        The epsilon conversion factor

    Returns
    ----------
    dihedral_types : list
        The dihedral types by number in the structure
    unique_dihedral_types : OrderedDict, ([(float, float, float, float, float, float, float, float,
                                            str, str, str, str, str, str, str, str), unique_number])
        This provides the unique dihedrals types, numbering, and the data values
        so it can easily be extracted.  An example of the OrderedDict is below:
        OrderedDict([(dihed_type_RB_c0, dihed_type_RB_c1, dihed_type_RB_c2,
                      dihed_type_RB_c3, dihed_type_RB_c4, dihed_type_RB_c5,
                      dihed_type_scee, dihed_type_scnb,
                      dihed_atom_1_type, dihed_atom_2_type,
                      dihed_atom_3_type, dihed_atom_4_type,
                      dihed_atom_1_res_type, dihed_atom_2_res_type,
                      dihed_atom_3_res_type, dihed_atom_4_res_type), 1,
                     ...,
                     (dihed_type_RB_c0, dihed_type_RB_c1, dihed_type_RB_c2,
                      dihed_type_RB_c3, dihed_type_RB_c4, dihed_type_RB_c5,
                      dihed_type_scee, dihed_type_scnb,
                      dihed_atom_1_type, dihed_atom_2_type,
                      dihed_atom_3_type, dihed_atom_4_type,
                      dihed_atom_1_res_type, dihed_atom_2_res_type,
                      dihed_atom_3_res_type, dihed_atom_4_res_type), n])
    """
    if use_rb_torsions:
        unique_dihedral_types = _get_unique_rb_torsion_types(
            structure, epsilon_conversion_factor
        )

        dihedral_types = [
            unique_dihedral_types[
                _get_dihedral_rb_torsion_key(
                    dihedral, epsilon_conversion_factor
                )
            ]
            for dihedral in structure.rb_torsions
        ]

    elif use_dihedrals:
        print_warn_text = (
            "WARNING: Using the charmm style and impropers is not "
            "available in the current version of this psf, pdb, and GOMC writer."
        )
        warn(print_warn_text)
        return None, None

    unique_dihedral_check_dict = OrderedDict()
    for i_value_dihed, i_key_dihed in unique_dihedral_types.items():
        i_value_duplicated = False
        for j_value_dihed, j_key_dihed in unique_dihedral_types.items():
            j_value_dihed_reorder = (
                j_value_dihed[0],
                j_value_dihed[1],
                j_value_dihed[2],
                j_value_dihed[3],
                j_value_dihed[4],
                j_value_dihed[5],
                j_value_dihed[6],
                j_value_dihed[7],
                j_value_dihed[11],
                j_value_dihed[10],
                j_value_dihed[9],
                j_value_dihed[8],
                j_value_dihed[15],
                j_value_dihed[14],
                j_value_dihed[13],
                j_value_dihed[12],
            )

            if i_value_dihed == j_value_dihed_reorder:
                i_value_duplicated = True
                if i_value_dihed[8] > j_value_dihed[8]:
                    unique_dihedral_check_dict.update(
                        {j_value_dihed: len(unique_dihedral_check_dict) + 1}
                    )
                else:
                    unique_dihedral_check_dict.update(
                        {i_value_dihed: len(unique_dihedral_check_dict) + 1}
                    )
        if i_value_duplicated is False:
            unique_dihedral_check_dict.update(
                {i_value_dihed: len(unique_dihedral_check_dict) + 1}
            )

    unique_dihedral_types = OrderedDict(
        [(y, x) for y, x in unique_dihedral_check_dict.items()]
    )

    return dihedral_types, unique_dihedral_types


def _get_impropers(structure, epsilon_conversion_factor):
    """
    Get a list of unique improper types that are used to create the Charmm style
    parameter file (i.e. force field file).  This also checks that the alternately
    ordered improper types are considered the same unique improper type.

    Parameters
    ----------
    structure : parmed.Structure
        The parmed structure
    epsilon_conversion_factor : float or int
        The epsilon conversion factor

    Returns
    ----------
    improper_types : list
        The improper types by number in the structure
    unique_improper_types : OrderedDict, ([(float, float, str, (str, str, str), str, (str, str, str)), unique_number])
        This provides the unique improper types, numbering, and the data values
        so it can easily be extracted.  An example of the OrderedDict is below:
        OrderedDict([(improper_k_constant, improper_psi_o,
          improper_atom_1_type, (improper_atom_2_type,
          improper_atom_3_type, improper_atom_4_type),
          improper_atom_1_res_type, (improper_atom_2_res_type,
          improper_atom_3_res_type, improper_atom_4_res_type)
         ), 1,
         ...,
         (improper_k_constant, improper_psi_o,
          improper_atom_1_type, (improper_atom_2_type,
          improper_atom_3_type, improper_atom_4_type),
          improper_atom_1_res_type, (improper_atom_2_res_type,
          improper_atom_3_res_type, improper_atom_4_res_type)
         ), n ])
    """
    unique_improper_types = _get_unique_improper_types(
        structure, epsilon_conversion_factor
    )
    improper_types = [
        unique_improper_types[
            _get_improper_type_key(improper, epsilon_conversion_factor)
        ]
        for improper in structure.impropers
    ]

    unique_improper_check_dict = OrderedDict()
    for i_value_improper, i_key_improper in unique_improper_types.items():
        i_value_duplicated = False

        i_value_improper_k_constant = i_value_improper[0]
        i_value_improper_psi_o = i_value_improper[1]

        i_value_impr_1_atoms_reorder = i_value_improper[2]
        i_value_impr_atoms_2_3_4_reorder = sorted(
            set(
                i_value_improper[3][0],
                i_value_improper[3][1],
                i_value_improper[3][2],
            )
        )

        i_value_impr_1_res_reorder = i_value_improper[4]
        i_value_impr_res_2_3_4_reorder = sorted(
            set(
                i_value_improper[5][0],
                i_value_improper[5][1],
                i_value_improper[5][2],
            )
        )
        i_improper_reformed = (
            i_value_improper_k_constant,
            i_value_improper_psi_o,
            i_value_impr_1_atoms_reorder,
            i_value_impr_atoms_2_3_4_reorder,
            i_value_impr_1_res_reorder,
            i_value_impr_res_2_3_4_reorder,
        )

        for j_value_improper, j_key_improper in unique_improper_types.items():
            j_value_improper_k_constant = j_value_improper[0]
            j_value_improper_psi_o = j_value_improper[1]

            j_value_impr_1_atoms_reorder = j_value_improper[2]
            j_value_impr_atoms_2_3_4_reorder = sorted(
                set(
                    j_value_improper[3][0],
                    j_value_improper[3][1],
                    j_value_improper[3][2],
                )
            )
            j_value_impr_1_res_reorder = j_value_improper[4]
            j_value_impr_res_2_3_4_reorder = sorted(
                set(
                    j_value_improper[5][0],
                    j_value_improper[5][1],
                    j_value_improper[5][2],
                )
            )

            j_improper_reformed = (
                j_value_improper_k_constant,
                j_value_improper_psi_o,
                j_value_impr_1_atoms_reorder,
                j_value_impr_atoms_2_3_4_reorder,
                j_value_impr_1_res_reorder,
                j_value_impr_res_2_3_4_reorder,
            )

            if i_improper_reformed == j_improper_reformed:
                i_value_duplicated = True
                if i_value_improper[2] > j_value_improper[2]:
                    unique_improper_check_dict.update(
                        {j_value_improper: len(unique_improper_check_dict) + 1}
                    )
                else:
                    unique_improper_check_dict.update(
                        {i_value_improper: len(unique_improper_check_dict) + 1}
                    )
        if i_value_duplicated is False:
            unique_improper_check_dict.update(
                {i_value_improper: len(unique_improper_check_dict) + 1}
            )

    unique_improper_types = OrderedDict(
        [(y, x) for y, x in unique_improper_check_dict.items()]
    )

    return improper_types, unique_improper_types


def unique_atom_naming(
    structure, residue_id_list, residue_names_list, bead_to_atom_name_dict=None
):
    """
    Generates unique atom/bead names for each molecule, which is required for some
    simulation types (Example: special Monte Carlo moves)

    Parameters
    ----------
    structure : compound object
    residue_id_list : list, in sequential order
            The residue ID for every atom in the system
    residue_names_list : list, in sequential order
        The atom names for every atom in the system
    bead_to_atom_name_dict: dictionary ; optional, default =None
        For all atom names/elements/beads with 2 or less digits, this converts
        the atom name in the GOMC psf and pdb files to a unique atom name,
        provided they do not exceed 3844 atoms (62^2) of the same name/element/bead
        per residue. For all atom names/elements/beads with 3 digits, this converts
        the atom name in the GOMC psf and pdb files to a unique atom name,
        provided they do not exceed 62 of the same name/element pre residue.
        Example dictionary: {'_CH3':'C', '_CH2':'C', '_CH':'C', '_HC':'C'}

    Returns
     ----------
    unique_individual_atom_names_dict : dictionary
        All the unique atom names comno_piled into a dictionary.
    individual_atom_names_list : list, in sequential  order
        The atom names for every atom in the system
    missing_bead_to_atom_name : list, in sequential  order
        The bead names of any atoms beads that did not have a name specificed to them
        via the bead_to_atom_name_dict
    """
    unique_individual_atom_names_dict = {}
    individual_atom_names_list = []
    missing_bead_to_atom_name = []
    for i, atom in enumerate(structure.atoms):
        interate_thru_names = True
        j = 0
        while interate_thru_names is True:
            j = j + 1
            if str(atom.name)[:1] == "_":
                if (
                    bead_to_atom_name_dict is not None
                    and (str(atom.name) in bead_to_atom_name_dict) is True
                ):
                    if len(bead_to_atom_name_dict[str(atom.name)]) > 2:
                        text_to_write = (
                            "ERROR: only enter atom names that have 2 or less digits"
                            + " in the Bead to atom naming dictionary (bead_to_atom_name_dict)."
                        )
                        warn(text_to_write)
                        return None, None, None
                    else:
                        atom_name_value = bead_to_atom_name_dict[str(atom.name)]
                        no_digits_atom_name = 2
                else:
                    missing_bead_to_atom_name.append(1)
                    atom_name_value = "BD"
                    no_digits_atom_name = 2
            elif len(str(atom.name)) > 2:
                if len(str(atom.name)) == 3:
                    no_digits_atom_name = 1
                    atom_name_value = atom.name
                else:
                    text_to_write = (
                        "ERROR: atom numbering will not work propery at"
                        + " the element has more than 4 charaters"
                    )
                    warn(text_to_write)
                    return None, None, None
            else:
                no_digits_atom_name = 2
                atom_name_value = atom.name
            atom_name_iteration = str(atom_name_value) + str(
                base10_to_base62_alph_num(j)
            )
            atom_res_no_resname_atomname_iteration = (
                str(residue_id_list[i])
                + "_"
                + str(residue_names_list[i])
                + "_"
                + atom_name_iteration
            )

            if (
                unique_individual_atom_names_dict.get(
                    str(atom_res_no_resname_atomname_iteration)
                )
                is None
            ):
                unique_individual_atom_names_dict.update(
                    {atom_res_no_resname_atomname_iteration: i + 1}
                )
                interate_thru_names = False
                individual_atom_names_list.append(
                    str(atom_name_value)
                    + str(
                        str(base10_to_base62_alph_num(j))[-no_digits_atom_name:]
                    )
                )

    if sum(missing_bead_to_atom_name) > 0:
        warn(
            "NOTE: All bead names were not found in the Bead to atom naming dictionary (bead_to_atom_name_dict) "
        )

    return [
        unique_individual_atom_names_dict,
        individual_atom_names_list,
        missing_bead_to_atom_name,
    ]


def _lengths_angles_to_vectors(lengths, angles, precision=6):
    """Converts the length and angles into CellBasisVectors

    Parameters
    ----------
    lengths : list-like, shape=(3,), dtype=float
        Lengths of the edges of the box (user chosen units).
    angles : list-like, shape=(3,), dtype=float, default=None
        Angles (in degrees) that define the tilt of the edges of the box. If
        None is given, angles are assumed to be [90.0, 90.0, 90.0]. These are
        also known as alpha, beta, gamma in the crystallography community.
    precision : int, optional, default=6
        Control the precision of the floating point representation of box
        attributes. If none provided, the default is 6 decimals.

    Returns
    -------
    box_vectors: numpy.ndarray, [[float, float, float], [float, float, float], [float, float, float]]
        Three (3) sets vectors for box 0 each with 3 float values, which represent
        the vectors for the Charmm-style systems (units are the same as entered for lengths)

    """

    (a, b, c) = lengths

    (alpha, beta, gamma) = np.deg2rad(angles)
    cos_a = np.clip(np.cos(alpha), -1.0, 1.0)
    cos_b = np.clip(np.cos(beta), -1.0, 1.0)
    cos_g = np.clip(np.cos(gamma), -1.0, 1.0)

    sin_a = np.clip(np.sin(alpha), -1.0, 1.0)
    sin_b = np.clip(np.sin(beta), -1.0, 1.0)
    sin_g = np.clip(np.sin(gamma), -1.0, 1.0)
    a_vec = np.asarray([a, 0.0, 0.0])

    b_x = b * cos_g
    b_y = b * sin_g
    b_vec = np.asarray([b_x, b_y, 0.0])

    c_x = c * cos_b
    c_cos_y_term = (cos_a - (cos_b * cos_g)) / sin_g
    c_y = c * c_cos_y_term
    c_z = c * np.sqrt(1 - np.square(cos_b) - np.square(c_cos_y_term))
    c_vec = np.asarray([c_x, c_y, c_z])
    box_vectors = np.asarray((a_vec, b_vec, c_vec))
    box_vectors.reshape(3, 3)
    # still leaves some floating values in some cases
    box_vectors = np.around(box_vectors, decimals=precision)

    return box_vectors


def _check_fixed_bonds_angles_lists(
    gomc_fix_bonds_and_or_angles,
    gomc_fix_bonds_and_or_angles_selection,
    residues,
):
    """Check the GOMC fixed bonds and angles lists for input errors.

    Parameters
    ----------
    gomc_fix_bonds_and_or_angles : list of strings, [str, ..., str]
        A list of the residues (i.e., molecules since GOMC currently considers a
        a whole molecule as a residue) to have their bonds and/or angles held
        rigid/fixed for the GOMC simulation engine.
        The `gomc_fix_bonds_angles`, `gomc_fix_bonds`, `gomc_fix_angles` are the only possible
        variables from the `Charmm` object to be entered.
        In GOMC, the residues currently are the same for every bead or atom in
        the molecules. Therefore, when the residue is selected, the whole molecule
        is selected.
    gomc_fix_bonds_and_or_angles_selection : str
        The name of the variable that is used but formatted as a string, which is fed
        to the error and information outputs. The
        `gomc_fix_bonds_angles`, `gomc_fix_bonds`, `gomc_fix_angles` are the only possible
        variables from the `Charmm` object to be entered.
        Whichever variable you choose, the variable name is just input as a
        string here. For example, if `gomc_fix_bonds_and_or_angles` is equal to
        gomc_fix_bonds_angles, then this should be 'gomc_fix_bonds_angles'
        (i.e., `gomc_fix_bonds_and_or_angles_selection` = 'gomc_fix_bonds_angles').
    residues : list, [str, ..., str]
        Labels of unique residues in the Compound. Residues are assigned by
        checking against Compound.name.  Only supply residue names as 4 character
        strings, as the residue names are truncated to 4 characters to fit in the
        psf and pdb file.

    Returns
    -------
    Provides a ValueError or TypeError if the input is not correct.
    """

    if gomc_fix_bonds_and_or_angles is not None and not isinstance(
        gomc_fix_bonds_and_or_angles, list
    ):
        print_error_message = (
            "ERROR: Please ensure the residue names in the ({}) variable "
            "are in a list.".format(gomc_fix_bonds_and_or_angles_selection)
        )
        raise TypeError(print_error_message)

    if isinstance(gomc_fix_bonds_and_or_angles, list):
        for gomc_fix_i in gomc_fix_bonds_and_or_angles:
            if gomc_fix_i not in residues:
                print_error_message = (
                    "ERROR: Please ensure that all the residue names in the "
                    "{} list are also in the residues list.".format(
                        gomc_fix_bonds_and_or_angles_selection
                    )
                )
                raise ValueError(print_error_message)
            elif not isinstance(gomc_fix_i, str):
                print_error_message = "ERROR: Please enter a fix_res_bonds list with only string values."
                raise TypeError(print_error_message)
            else:
                print(
                    "INFORMATION: The following residues will have these fixed parameters: "
                    + "gomc_fix_bonds = {}".format(gomc_fix_bonds_and_or_angles)
                )


# Currently the NBFIX is disabled as since only the OPLS and TRAPPE force fields are currently supported
class Charmm:
    """Generates a Charmm object that is required to produce the Charmm style parameter
    (force field), PDB, PSF files, which are usable in the GOMC and NAMD engines.
    Additionally, this Charmm object is also used in generating the GOMC control file.

    The units for the GOMC data files.
        * Mw = g/mol
        * charge = e
        * Harmonic bonds : Kb = kcal/mol, b0 = Angstroms
        * Harmonic angles : Ktheta = kcal/mole/rad**2 , Theta0 = degrees
        * Dihedral angles: Ktheta = kcal/mole, n = interger (unitless), delta = degrees
        * Improper angles (currently unavailable) : TBD
        * LJ-NONBONDED : epsilon = kcal/mol, Rmin/2 = Angstroms
        * Mie-NONBONDED (currently unavailable): epsilon = K, sigma = Angstroms, n = interger (unitless)
        * Buckingham-NONBONDED (currently unavailable): epsilon = K, sigma = Angstroms, n = interger (unitless)
        * LJ-NBFIX (currently unavailable) : epsilon = kcal/mol, Rmin = Angstroms
        * Mie-NBFIX (currently unavailable) : same as Mie-NONBONDED
        * Buckingham-NBFIX (currently unavailable) : same as Buckingham-NONBONDED

    Note: units are the same as the NAMD units and the LAMMPS real units.  The atom style
    is the same as the lammps 'full' atom style format.

    Parameters
    ----------
    structure_box_0 : mbuild Compound object (mbuild.Compound) or mbuild Box object (mbuild.Box);
        If the structure has atoms/beads it must be an mbuild Compound.
        If the structure is empty it must be and mbuild Box object.
        Note: If 1 structures are provided (i.e., only structure_box_0),
        it must be an mbuild Compound.
        Note: If 2 structures are provided,
        only 1 structure can be an empty box (i.e., either structure_box_0 or structure_box_1)
    filename_box_0 : str
        The file name of the output file for structure_box_0.  Note: the extension should
        not be provided, as multiple extension (.pdb and .psf) are added to this name.
    structure_box_1 : mbuild Compound object (mbuild.Compound) or mbuild Box object (mbuild.Box), default = None;
        If the structure has atoms/beads it must be an mbuild Compound.
        Note: When running a GEMC or GCMC simulation the box 1 stucture should be input
        here.  Otherwise, there is no guarantee that any of the atom type and force field
        information will all work together correctly with box 0, if it is built separately.
        Note: If 2 structures are provided, only 1 structure can be an empty box
        (i.e., either structure_box_0 or structure_box_1).
    filename_box_1 : str , default = None
        The file name of the output file for structure_box_1 (Ex: for GCMC or GEMC simulations
        which have mulitiple simulation boxes).  Note: the extension should
        not be provided, as multiple extension (.pdb and .psf) are added to this name.
        Note: When running a GEMC or GCMC simulation the box 1 stucture should be input
        here.  Otherwise, there is no guarantee that any of the atom type and force field
        information will all work together correctly with box 0, if it is built separately.
    non_bonded_type : str, default = 'LJ' (i.e., Lennard-Jones )
        Specify the type of non-bonded potential for the GOMC force field files.
        Note: Currently, on the 'LJ' potential is supported.
    residues : list, [str, ..., str]
        Labels of unique residues in the Compound. Residues are assigned by
        checking against Compound.name.  Only supply residue names as 4 character
        strings, as the residue names are truncated to 4 characters to fit in the
        psf and pdb file.
    forcefield_selection : str or dictionary, default = None
        Apply a forcefield to the output file by selecting a force field XML file with
        its path or by using the standard force field name provided the `foyer` package.
        Note: to write the NAMD/GOMC force field, pdb, and psf files, the
        residues and forcefields must be provided in a str or
        dictionary.  If a dictionary is provided all residues must
        be specified to a force field.
        * Example dict for FF file: {'ETH' : 'oplsaa.xml', 'OCT': 'path_to_file/trappe-ua.xml'}

        * Example str for FF file: 'path_to file/trappe-ua.xml'

        * Example dict for standard FF names : {'ETH' : 'oplsaa', 'OCT': 'trappe-ua'}

        * Example str for standard FF names: 'trappe-ua'

        * Example of a mixed dict with both : {'ETH' : 'oplsaa', 'OCT': 'path_to_file/'trappe-ua.xml'}
    detect_forcefield_style : boolean, default = True
        If True, format NAMD/GOMC/LAMMPS parameters based on the contents of
        the parmed structure_box_0 and structure_box_1.
    gomc_fix_bonds_angles : list, default = None
        When list of residues is provided, the selected residues will have
        their bonds and angles fixed in the GOMC engine.  This is specifically
        for the GOMC engine and it changes the residue's bond constants (Kbs)
        and angle constants (Kthetas) values to 999999999999 in the
        FF file (i.e., the .inp file).
    bead_to_atom_name_dict : dict, optional, default =None
        For all atom names/elements/beads with 2 or less digits, this converts
        the atom name in the GOMC psf and pdb files to a unique atom name,
        provided they do not exceed 3844 atoms (62^2) of the same name/element/bead
        per residue. For all atom names/elements/beads with 3 digits, this converts
        the atom name in the GOMC psf and pdb files to a unique atom name,
        provided they do not exceed 62 of the same name/element pre residue.

        * Example dictionary: {'_CH3':'C', '_CH2':'C', '_CH':'C', '_HC':'C'}

        * Example name structure: {atom_type: first_part_pf atom name_without_numbering}

    fix_residue : list  or None, default = None
        Changes occcur in the pdb file only.
        When residues are listed here, all the atoms in the residue are
        fixed and can not move via setting the Beta values in the PDB
        file to 1.00.
        If neither fix_residue or fix_residue_in_box lists a
        residue or both equal None, then the Beta values for all the atoms
        in the residue are free to move in the simulation and Beta values
        in the PDB file is set to 0.00
    fix_residue_in_box : list  or None, default = None
        Changes occcur in the pdb file only.
        When residues are listed here, all the atoms in the residue
        can move within the box but cannot be transferred between boxes
        via setting the Beta values in the PDB file to 2.00.
        If neither fix_residue or fix_residue_in_box lists a
        residue or both equal None, then the Beta values for all the atoms
        in the residue are free to move in the simulation and Beta values
        in the PDB file is set to 0.00
    ff_filename : str, default =None
        If a string, it will write the  force field files that work in
        GOMC and NAMD structures.
    reorder_res_in_pdb_psf : bool, default =False
        If False, the order of of the atoms in the pdb file is kept in
        its original order, as in the Compound sent to the writer.
        If True, the order of the atoms is reordered based on their
        residue names in the 'residues' list that was entered.

    Attributes
    ----------
    input_error : bool
        This error is typically incurred from an error in the user's input values.
        However, it could also be due to a bug, provided the user is inputting
        the data as this Class intends.
    structure_box_0 : mbuild.compound.Compound
        The mbuild Compound for the input box 0
    structure_box_1 : mbuild.compound.Compound or None, default = None
        The mbuild Compound for the input box 1
    filename_box_0 : str
        The file name of the output file for structure_box_0.  Note: the extension should
        not be provided, as multiple extension (.pdb and .psf) are added to this name.
    filename_box_1 : str or None , default = None
        The file name of the output file for structure_box_1.  Note: the extension should
        not be provided, as multiple extension (.pdb and .psf) are added to this name.
        (i.e., either structure_box_0 or structure_box_1).
    non_bonded_type : str, default = 'LJ' (i.e., Lennard-Jones )
        Specify the type of non-bonded potential for the GOMC force field files.
        Note: Currently, on the 'LJ' potential is supported.
    residues : list, [str, ..., str]
        Labels of unique residues in the Compound. Residues are assigned by
        checking against Compound.name.  Only supply residue names as 4 character
        strings, as the residue names are truncated to 4 characters to fit in the
        psf and pdb file.
    forcefield_selection : str or dictionary, default = None
        Apply a forcefield to the output file by selecting a force field XML file with
        its path or by using the standard force field name provided the `foyer` package.
        Note: to write the NAMD/GOMC force field, pdb, and psf files, the
        residues and forcefields must be provided in a str or
        dictionary.  If a dictionary is provided all residues must
        be specified to a force field.

        * Example dict for FF file: {'ETH' : 'oplsaa.xml', 'OCT': 'path_to_file/trappe-ua.xml'}

        * Example str for FF file: 'path_to file/trappe-ua.xml'

        * Example dict for standard FF names : {'ETH' : 'oplsaa', 'OCT': 'trappe-ua'}

        * Example str for standard FF names: 'trappe-ua'

        * Example of a mixed dict with both : {'ETH' : 'oplsaa', 'OCT': 'path_to_file/'trappe-ua.xml'}

    detect_forcefield_style : bool, default = True
        If True, format NAMD/GOMC/LAMMPS parameters based on the contents of
        the parmed structure_box_0 and structure_box_1
    gomc_fix_bonds_angles : list, default = None
        When list of residues is provided, the selected residues will have
        their bonds and angles fixed and will ignore the relative bond energies and
        related angle energies in the GOMC engine. Note that GOMC
        does not sample bond stretching. This is specifically
        for the GOMC engine and it changes the residue's bond constants (Kbs)
        and angle constants (Kthetas) values to 999999999999 in the
        FF file (i.e., the .inp file).
        If the residues are listed in either the gomc_fix_angles or the gomc_fix_bonds_angles
        lists, the angles will be fixed for that residue.
        If the residues are listed in either the gomc_fix_bonds or the gomc_fix_bonds_angles
        lists, the bonds will be fixed for that residue.
        NOTE if this option is utilized it may cause issues if using the FF file in NAMD.
    gomc_fix_bonds : list, default = None
        When list of residues is provided, the selected residues will have their
        relative bond energies ignored in the GOMC engine. Note that GOMC
        does not sample bond stretching. This is specifically
        for the GOMC engine and it changes the residue's bond constants (Kbs)
        values to 999999999999 in the FF file (i.e., the .inp file).
        If the residues are listed in either the gomc_fix_bonds or the gomc_fix_bonds_angles
        lists, the relative bond energy will be ignored.
        NOTE if this option is utilized it may cause issues if using the FF file in NAMD.
    gomc_fix_angles : list, default = None
        When list of residues is provided, the selected residues will have
        their angles fixed and will ignore the related angle energies in the GOMC engine.
        This is specifically for the GOMC engine and it changes the residue's angle
        constants (Kthetas) values to 999999999999 in the FF file (i.e., the .inp file),
        which fixes the angles and ignores related angle energy.
        If the residues are listed in either the gomc_fix_angles or the gomc_fix_bonds_angles
        lists, the angles will be fixed and the related angle energy will be ignored
        for that residue.
        NOTE if this option is utilized it may cause issues if using the FF file in NAMD.
    bead_to_atom_name_dict : dict, optional, default =None
        For all atom names/elements/beads with 2 or less digits, this converts
        the atom name in the GOMC psf and pdb files to a unique atom name,
        provided they do not exceed 3844 atoms (62^2) of the same name/element/bead
        per residue. For all atom names/elements/beads with 3 digits, this converts
        the atom name in the GOMC psf and pdb files to a unique atom name,
        provided they do not exceed 62 of the same name/element pre residue.

        * Example dictionary: {'_CH3':'C', '_CH2':'C', '_CH':'C', '_HC':'C'}

        * Example name structure: {atom_type: first_part_pf atom name_without_numbering}

    fix_residue : list  or None, default = None
        Changes occcur in the pdb file only.
        When residues are listed here, all the atoms in the residue are
        fixed and can not move via setting the Beta values in the PDB
        file to 1.00.
        If neither fix_residue or fix_residue_in_box lists a
        residue or both equal None, then the Beta values for all the atoms
        in the residue are free to move in the simulation and Beta values
        in the PDB file is set to 0.00
    fix_residue_in_box : list  or None, default = None
        Changes occcur in the pdb file only.
        When residues are listed here, all the atoms in the residue
        can move within the box but cannot be transferred between boxes
        via setting the Beta values in the PDB file to 2.00.
        If neither fix_residue or fix_residue_in_box lists a
        residue or both equal None, then the Beta values for all the atoms
        in the residue are free to move in the simulation and Beta values
        in the PDB file is set to 0.00
    ff_filename : str, default =None
        If a string, it will write the  force field files that work in
        GOMC and NAMD structures.
    reorder_res_in_pdb_psf : bool, default =False
        If False, the order of of the atoms in the pdb file is kept in
        its original order, as in the Compound sent to the writer.
        If True, the order of the atoms is reordered based on their
        residue names in the 'residues' list that was entered.
    box_0 : Box
        The Box class that contains the attributes Lx, Ly, Lz for the length
        of the box 0 (units in nanometers (nm)). It also contains the xy, xz, and yz Tilt factors
        needed to displace an orthogonal box's xy face to its
        parallelepiped structure for box 0.
    box_1 : Box
        The Box class that contains the attributes Lx, Ly, Lz for the length
        of the box 1 (units in nanometers (nm)). It also contains the xy, xz, and yz Tilt factors
        needed to displace an orthogonal box's xy face to its
        parallelepiped structure for box 0.
    box_0_vectors : numpy.ndarray, [[float, float, float], [float, float, float], [float, float, float]]
        Three (3) sets vectors for box 0 each with 3 float values, which represent
        the vectors for the Charmm-style systems (units in Angstroms (Ang))
    box_1_vectors : numpy.ndarray, [[float, float, float], [float, float, float], [float, float, float]]
        Three (3) sets vectors for box 1 each with 3 float values, which represent
        the vectors for the Charmm-style systems (units in Angstroms (Ang))
    structure_box_0_ff : parmed.structure.Structure
        The box 0 structure (structure_box_0) after all the provided
        force fields are applied.
    structure_box_1_ff : parmed.structure.Structure
        The box 0 structure (structure_box_0) after all the provided
        force fields are applied. This only exists if the box 1 structure
        (structure_box_1) is provided.
    coulomb14scalar_dict_box_0 : dict
        The residue/moleclues (key) of box 0 and their corresponding
        coulombic 1-4 scalers (value).  Note: NAMD and GOMC can only have one (1)
        value for the coulombic 1-4 scalers, as they both only accept a
        single value in the NAMD and GOMC control files.
    coulomb14scalar_dict_box_1 : dict
        The residue/moleclues  (key) of box 1 and their corresponding
        coulombic 1-4 scalers (value).  Note: NAMD and GOMC can only have one (1)
        value for the coulombic 1-4 scalers, as they both only accept a
        single value in the NAMD and GOMC control files.
        This only exists if the box 1 structure (structure_box_1) is provided.
    LJ14scalar_dict_box_0 : dict
        The residue/moleclues (key) of box 0 and their corresponding
        Lennard-Jones (LJ) 1-4 scalers (value).  Note: NAMD and GOMC can have
        multiple values for the LJ 1-4 scalers, since they are provided as an
        individual input for each atom type in the force field (.inp) file.
    LJ14scalar_dict_box_1 : dict
        The residue/moleclues (key) of box 1 and their corresponding
        Lennard-Jones (LJ) 1-4 scalers (value).  Note: NAMD and GOMC can have
        multiple values for the LJ 1-4 scalers, since they are provided as an
        individual input for each atom type in the force field (.inp) file.
        This only exists if the box 1 structure (structure_box_1) is provided.
    residues_applied_list_box_0 : list
        The residues in box 0 that were found and had the force fields applied to them.
    residues_applied_list_box_1 : list
        The residues in box 1 that were found and had the force fields applied to them.
        This only exists if the box 1 structure (structure_box_1) is provided.
    boxes_for_simulation : int, [0, 1]
        The number of boxes used when writing the Charmm object and force fielding
        the system. If only box 0 is provided, the value is 0. If box 0 and box 1
        are provided, the value is 1.
    epsilon_dict : dict {str: float or int}
        The uniquely numbered atom type (key) and it's non-bonded epsilon
        coefficient (value).
    sigma_dict : dict {str: float or int}
        The uniquely numbered atom type (key) and it's non-bonded sigma
        coefficient (value).
    LJ_1_4_dict : dict {str: float or int}
        The uniquely numbered atom type (key) and it's non-bonded 1-4
        Lennard-Jones, LJ, scaling factor (value).
    coul_1_4 : float or int
        The non-bonded 1-4 coulombic scaling factor, which is the
        same for all the residues/molecules, regardless if
        differenct force fields are utilized.  Note: if
        1-4 coulombic scaling factor is not the same for all
        molecules the Charmm object will fail with an error.
    combined_1_4_coul_dict_per_residue : dict, {str: float or int}
        The residue name/molecule (key) and it's non-bonded 1-4 coulombic
        scaling factor (value).
    forcefield_dict : dict
        The uniquely numbered atom type (key) with it's corresponding
        foyer atom typing and residue name.  The residue name is added
        to provide a distinction between other residues with the same
        atom types.  This allows the CHARMM force field to fix the
        bonds and angles specific residues without effecting other
        residues with the same atom types.
    all_individual_atom_names_list : list
        A list of all the atom names for the combined structures
        (box 0 and box 1 (if supplied)), in order.
    all_residue_names_List : list
        A list of all the residue names for the combined structures
        (box 0 and box 1 (if supplied)), in order.
    max_residue_no : int
        The maximum number that the residue number will count to
        before restarting the counting back to 1, which is predetermined
        by the PDB format. This is a constant, which equals 9999
    max_resname_char : int
        The maximum number of characters allowed in the residue name,
        which is predetermined by the PDB format. This is a constant,
        which equals 4.
    all_res_unique_atom_name_dict : dict, {str : [str, ..., str]}
        A dictionary that provides the residue names (keys) and a list
        of the unique atom names in the residue (value), for the
        combined structures (box 0 and box 1 (if supplied)).

    Notes
    -----
    Impropers, Urey-Bradleys, and NBFIX are not currenly supported.
    Currently the NBFIX is disabled as since only the OPLS and TRAPPE force fields are supported.
    OPLS and CHARMM forcefield styles are supported (without impropers),
    AMBER forcefield styles are NOT supported.

    The atom typing is currently provided via a base 52 numbering (capital and lowercase lettering).
    This base 52 numbering allows for (52)^4 unique atom types.

    Unique atom names are provided if the system do not exceed 3844 atoms (62^2) of the same
    name/bead per residue (base 62 numbering). For all atom names/elements with 3 or less digits,
    this converts the atom name in the GOMC psf and pdb files to a unique atom name,
    provided they do not exceed 62 of the same name/element pre residue.

    Generating an empty box (i.e., pdb and psf files):
    Single Box system: Enter residues = [], but the accompanying structure (structure_box_0)
    must be an empty mb.Box. However, when doing this, the forcefield_selection
    must be supplied, or it will provide an error
    (i.e., forcefield_selection can not be equal to None).
    Dual Box System: Enter an empty mb.Box structure for either structure_box_0 or
    structure_box_1.

    In this current FF/psf/pdb writer, a residue type is essentially a molecule type.
    Therefore, it can only correctly write systems where every bead/atom in the molecule
    has the same residue name, and the residue name is specific to that molecule type.
    For example: a protein molecule with many residue names is not currently supported,
    but is planned to be supported in the future.
    """

    def __init__(
        self,
        structure_box_0,
        filename_box_0,
        structure_box_1=None,
        filename_box_1=None,
        non_bonded_type="LJ",
        forcefield_selection=None,
        residues=None,
        detect_forcefield_style=True,
        gomc_fix_bonds_angles=None,
        gomc_fix_bonds=None,
        gomc_fix_angles=None,
        bead_to_atom_name_dict=None,
        fix_residue=None,
        fix_residue_in_box=None,
        ff_filename=None,
        reorder_res_in_pdb_psf=False,
    ):

        # set all input variables to the class
        self.structure_box_0 = structure_box_0
        self.filename_box_0 = filename_box_0
        self.structure_box_1 = structure_box_1
        self.filename_box_1 = filename_box_1
        self.non_bonded_type = non_bonded_type
        self.forcefield_selection = forcefield_selection
        self.residues = residues
        self.detect_forcefield_style = detect_forcefield_style
        self.gomc_fix_bonds_angles = gomc_fix_bonds_angles
        self.gomc_fix_bonds = gomc_fix_bonds
        self.gomc_fix_angles = gomc_fix_angles
        self.bead_to_atom_name_dict = bead_to_atom_name_dict
        self.fix_residue = fix_residue
        self.fix_residue_in_box = fix_residue_in_box
        self.ff_filename = ff_filename
        self.reorder_res_in_pdb_psf = reorder_res_in_pdb_psf

        # value to check for errors, with  self.input_error = True or False. Set to False initally
        self.input_error = False

        if not isinstance(self.structure_box_0, (Compound, Box)):
            self.input_error = True
            print_error_message = (
                "ERROR: The structure_box_0 expected to be of type: "
                "{} or {}, received: {}".format(
                    type(Compound()),
                    type(Box(lengths=[1, 1, 1])),
                    type(structure_box_0),
                )
            )
            raise TypeError(print_error_message)

        if self.structure_box_1 is not None and not isinstance(
            self.structure_box_1, (Compound, Box)
        ):
            self.input_error = True
            print_error_message = (
                "ERROR: The structure_box_1 expected to be of type: "
                "{} or {}, received: {}".format(
                    type(Compound()),
                    type(Box(lengths=[1, 1, 1])),
                    type(structure_box_1),
                )
            )
            raise TypeError(print_error_message)

        if isinstance(self.structure_box_0, Box) and isinstance(
            self.structure_box_1, Box
        ):
            self.input_error = True
            print_error_message = (
                "ERROR: Both structure_box_0 and structure_box_0 are empty Boxes {}. "
                "At least 1 structure must be an mbuild compound {} with 1 "
                "or more atoms in it".format(
                    type(Box(lengths=[1, 1, 1])), type(Compound())
                )
            )
            raise TypeError(print_error_message)

        if self.structure_box_1 is None and not isinstance(
            self.structure_box_0, Compound
        ):
            self.input_error = True
            print_error_message = (
                "ERROR: Only 1 structure is provided and it can not be an empty mbuild Box {}. "
                "it must be an mbuild compound {} with at least 1 "
                "or more atoms in it.".format(
                    type(Box(lengths=[1, 1, 1])), type(Compound())
                )
            )
            raise TypeError(print_error_message)

        if not isinstance(self.residues, list):
            self.input_error = True
            print_error_message = "ERROR: Please enter the residues list (residues) in a list format."
            raise TypeError(print_error_message)

        if isinstance(self.residues, list):
            for each_residue in self.residues:
                if not isinstance(each_residue, str):
                    self.input_error = True
                    print_error_message = "ERROR: Please enter a residues list (residues) with only string values."
                    raise TypeError(print_error_message)

        if self.residues is None:
            self.input_error = True
            print_error_message = (
                "ERROR: Please enter the residues list (residues)"
            )
            raise TypeError(print_error_message)
        if not isinstance(self.filename_box_0, str):
            self.input_error = True
            print_error_message = (
                "ERROR: Please enter the filename_box_0 as a string."
            )
            raise TypeError(print_error_message)

        unique_residue_test_name_list = []
        for res_m in range(0, len(self.residues)):
            if self.residues[res_m] not in unique_residue_test_name_list:
                unique_residue_test_name_list.append(self.residues[res_m])
        if len(unique_residue_test_name_list) != len(self.residues):
            self.input_error = True
            print_error_message = "ERROR: Please enter the residues list (residues) that has only unique residue names."
            raise ValueError(print_error_message)

        if self.filename_box_1 is not None and not isinstance(
            self.filename_box_1, str
        ):
            self.input_error = True
            print_error_message = (
                "ERROR: Please enter the filename_box_1 as a string."
            )
            raise TypeError(print_error_message)

        if self.ff_filename is not None:
            if not isinstance(self.ff_filename, str):
                self.input_error = True
                print_error_message = "ERROR: Please enter GOMC force field name (ff_filename) as a string."
                raise TypeError(print_error_message)
            if isinstance(self.ff_filename, str):
                extension_ff_name = os.path.splitext(self.ff_filename)[-1]
                if extension_ff_name == "":
                    self.ff_filename = self.ff_filename + ".inp"
                elif extension_ff_name == ".inp":
                    self.ff_filename = self.ff_filename + ""
                elif extension_ff_name != ".inp":
                    self.input_error = True
                    print_error_message = (
                        "ERROR: Please enter GOMC force field name without an "
                        "extention or the .inp extension."
                    )
                    raise ValueError(print_error_message)

        if self.forcefield_selection is not None:
            print(
                "write_gomcdata: forcefield_selection = "
                + str(self.forcefield_selection)
                + ", "
                + "residues = "
                + str(self.residues)
            )
            if not isinstance(self.forcefield_selection, (dict, str)):
                self.input_error = True
                print_error_message = (
                    "ERROR: The force field selection (forcefield_selection) "
                    "is not a string or a dictionary with all the residues specified "
                    'to a force field. -> String Ex: "path/trappe-ua.xml" or Ex: "trappe-ua" '
                    "Otherise provided a dictionary with all the residues specified "
                    "to a force field "
                    '->Dictionary Ex: {"Water" : "oplsaa", "OCT": "path/trappe-ua.xml"}, '
                    "Note: the file path must be specified the force field file if "
                    "a standard foyer force field is not used."
                )
                raise TypeError(print_error_message)

            if isinstance(self.forcefield_selection, str):
                ff_name = self.forcefield_selection
                self.forcefield_selection = {}
                for i in range(0, len(self.residues)):
                    self.forcefield_selection.update(
                        {self.residues[i]: ff_name}
                    )
                print(
                    "FF forcefield_selection = "
                    + str(self.forcefield_selection)
                )

        elif self.forcefield_selection is None:
            self.input_error = True
            print_error_message = "ERROR: Please enter the forcefield_selection as it was not provided."
            raise TypeError(print_error_message)

        if self.residues is not None and not isinstance(self.residues, list):
            self.input_error = True
            print_error_message = (
                "ERROR:  Please enter the residues (residues) in a list format"
            )
            raise TypeError(print_error_message)

        _check_fixed_bonds_angles_lists(
            self.gomc_fix_bonds_angles, "gomc_fix_bonds_angles", self.residues
        )
        _check_fixed_bonds_angles_lists(
            self.gomc_fix_bonds, "gomc_fix_bonds", self.residues
        )
        _check_fixed_bonds_angles_lists(
            self.gomc_fix_angles, "gomc_fix_angles", self.residues
        )

        if self.fix_residue is not None and not isinstance(
            self.fix_residue, list
        ):
            self.input_error = True
            print_error_message = (
                "ERROR: Please enter the fix_residue in a list format"
            )
            raise TypeError(print_error_message)

        if isinstance(self.fix_residue, list):
            for fix_residue_q in self.fix_residue:
                if fix_residue_q not in self.residues:
                    self.input_error = True
                    print_error_message = (
                        "Error: Please ensure that all the residue names in the fix_residue "
                        "list are also in the residues list."
                    )
                    raise ValueError(print_error_message)

        if self.fix_residue_in_box is not None and not isinstance(
            self.fix_residue_in_box, list
        ):
            self.input_error = True
            print_error_message = (
                "ERROR: Please enter the fix_residue_in_box in a list format."
            )
            raise TypeError(print_error_message)

        if isinstance(self.fix_residue_in_box, list):
            for fix_residue_in_box_q in self.fix_residue_in_box:
                if fix_residue_in_box_q not in self.residues:
                    self.input_error = True
                    print_error_message = (
                        "Error: Please ensure that all the residue names in the "
                        "fix_residue_in_box list are also in the residues list."
                    )
                    raise ValueError(print_error_message)

        if self.bead_to_atom_name_dict is not None and not isinstance(
            self.bead_to_atom_name_dict, dict
        ):
            self.input_error = True
            print_error_message = (
                "ERROR: Please enter the a bead type to atom in the dictionary "
                "(bead_to_atom_name_dict) so GOMC can properly evaluate the unique atom names"
            )
            raise TypeError(print_error_message)

        if isinstance(self.bead_to_atom_name_dict, dict):
            dict_list = []
            for key in self.bead_to_atom_name_dict.keys():
                dict_list.append(key)

            for dict_lis_i in dict_list:
                if not isinstance(dict_lis_i, str) or not isinstance(
                    self.bead_to_atom_name_dict[dict_lis_i], str
                ):
                    print_error_message = "ERROR: Please enter the bead_to_atom_name_dict with only string inputs."
                    raise TypeError(print_error_message)

        print("******************************")
        print("")

        self.sub_1_for_base_52 = 1

        if self.structure_box_1:
            self.boxes_for_simulation = 2
        else:
            self.boxes_for_simulation = 1

        # write the Force fields
        self.combined_1_4_lj_dict_per_residue = {}
        self.combined_1_4_coul_dict_per_residue = {}

        if self.structure_box_1:

            print(
                "GOMC FF writing each residues FF as a group for structure_box_0"
            )
            [
                self.structure_box_0_ff,
                self.coulomb14scalar_dict_box_0,
                self.LJ14scalar_dict_box_0,
                self.residues_applied_list_box_0,
            ] = specific_ff_to_residue(
                self.structure_box_0,
                forcefield_selection=self.forcefield_selection,
                residues=self.residues,
                reorder_res_in_pdb_psf=self.reorder_res_in_pdb_psf,
                boxes_for_simulation=self.boxes_for_simulation,
            )

            print(
                "GOMC FF writing each residues FF as a group for  structure_box_1"
            )
            [
                self.structure_box_1_ff,
                self.coulomb14scalar_dict_box_1,
                self.LJ14scalar_dict_box_1,
                self.residues_applied_list_box_1,
            ] = specific_ff_to_residue(
                self.structure_box_1,
                forcefield_selection=self.forcefield_selection,
                residues=self.residues,
                reorder_res_in_pdb_psf=self.reorder_res_in_pdb_psf,
                boxes_for_simulation=self.boxes_for_simulation,
            )

            self.structure_box_0_and_1_ff = (
                self.structure_box_0_ff + self.structure_box_1_ff
            )
            self.combined_1_4_lj_dict_per_residue.update(
                self.LJ14scalar_dict_box_0
            )
            self.combined_1_4_lj_dict_per_residue.update(
                self.LJ14scalar_dict_box_1
            )
            self.combined_1_4_coul_dict_per_residue.update(
                self.coulomb14scalar_dict_box_0
            )
            self.combined_1_4_coul_dict_per_residue.update(
                self.coulomb14scalar_dict_box_1
            )

            self.residues_applied_list_box_0_and_1 = (
                self.residues_applied_list_box_0
            )
            for res_iter in self.residues_applied_list_box_1:
                if res_iter not in self.residues_applied_list_box_0:
                    self.residues_applied_list_box_0_and_1.append(res_iter)

            for res_iter_0_1 in self.residues_applied_list_box_0_and_1:
                if res_iter_0_1 not in self.residues:
                    self.input_error = True
                    print_error_message = "ERROR: All the residues were not used from the forcefield_selection "
                    "string or dictionary.  There may be residues below other specified "
                    "residues in the mbuild.Compound hierarchy.  If so, the residues "
                    "acquire the residue's force fields, which is at the top of the "
                    "hierarchy.  Alternatively, residues that are not in the structure "
                    r"may have been specified."
                    raise ValueError(print_error_message)

            for res_iter_0_1 in self.residues:
                if res_iter_0_1 not in self.residues_applied_list_box_0_and_1:
                    self.input_error = True
                    print_error_message = (
                        "ERROR: All the residues were not used from the forcefield_selection "
                        "string or dictionary.  There may be residues below other specified "
                        "residues in the mbuild.Compound hierarchy.  If so, the residues "
                        "acquire the residue's force fields, which is at the top of the "
                        "hierarchy.  Alternatively, residues that are not in the structure "
                        "may have been specified."
                    )
                    raise ValueError(print_error_message)

            total_charge = sum(
                [atom.charge for atom in self.structure_box_0_ff]
            )
            if round(total_charge, 4) != 0.0:
                warn(
                    "System is not charge neutral for structure_box_0. "
                    "Total charge is {}.".format(total_charge)
                )

            total_charge = sum(
                [atom.charge for atom in self.structure_box_1_ff]
            )
            if round(total_charge, 4) != 0.0:
                warn(
                    "System is not charge neutral for structure_box_1. "
                    "Total charge is {}.".format(total_charge)
                )

            total_charge = sum(
                [atom.charge for atom in self.structure_box_0_and_1_ff]
            )
            if round(total_charge, 4) != 0.0:
                warn(
                    "System is not charge neutral for structure_0_and_1. "
                    "Total charge is {}.".format(total_charge)
                )

        else:
            print(
                "GOMC FF writing each residues FF as a group for structure_box_0"
            )
            [
                self.structure_box_0_ff,
                self.coulomb14scalar_dict_box_0,
                self.LJ14scalar_dict_box_0,
                self.residues_applied_list_box_0,
            ] = specific_ff_to_residue(
                self.structure_box_0,
                forcefield_selection=self.forcefield_selection,
                residues=self.residues,
                reorder_res_in_pdb_psf=self.reorder_res_in_pdb_psf,
                boxes_for_simulation=self.boxes_for_simulation,
            )

            self.combined_1_4_lj_dict_per_residue.update(
                self.LJ14scalar_dict_box_0
            )
            self.combined_1_4_coul_dict_per_residue.update(
                self.coulomb14scalar_dict_box_0
            )

            for res_iter_0 in self.residues_applied_list_box_0:
                if res_iter_0 not in self.residues:
                    self.input_error = True
                    print_error_message = (
                        "ERROR: All the residues were not used from the forcefield_selection "
                        "string or dictionary.  There may be residues below other specified "
                        "residues in the mbuild.Compound hierarchy.  If so, the residues "
                        "acquire the residue's force fields, which is at the top of the "
                        "hierarchy.  Alternatively, residues that are not in the structure "
                        "may have been specified."
                    )
                    raise ValueError(print_error_message)

            for res_iter_0 in self.residues:
                if res_iter_0 not in self.residues_applied_list_box_0:
                    self.input_error = True
                    print_error_message = (
                        "ERROR: All the residues were not used from the forcefield_selection "
                        "string or dictionary.  There may be residues below other specified "
                        "residues in the mbuild.Compound hierarchy.  If so, the residues "
                        "acquire the residue's force fields, which is at the top of the "
                        "hierarchy.  Alternatively, residues that are not in the structure "
                        "may have been specified."
                    )
                    raise ValueError(print_error_message)

            total_charge = sum(
                [atom.charge for atom in self.structure_box_0_ff]
            )
            if round(total_charge, 4) != 0.0:
                warn(
                    "System is not charge neutral for structure_box_0. "
                    "Total charge is {}.".format(total_charge)
                )

        print(
            "forcefield type from compound = " + str(self.forcefield_selection)
        )
        print(
            "coulomb14scale from compound = "
            + str(self.combined_1_4_coul_dict_per_residue)
        )
        print(
            "lj14scale from compound = "
            + str(self.combined_1_4_lj_dict_per_residue)
        )

        # lock the atom_style and unit_style for GOMC. Can be inserted into variables
        # once more functionality is built in
        self.atom_style = "full"
        self.unit_style = "real"
        # functional form type default.  Can be inserted into variables once more functionality is built in
        use_rb_torsions = True
        use_dihedrals = False
        use_urey_bradleys = False

        # Convert coordinates to real or other units (real only current option)
        if self.unit_style == "real":
            self.sigma_conversion_factor = 1
            self.epsilon_conversion_factor = 1
            self.mass_conversion_factor = 1

        if self.structure_box_1:
            self.types = np.array(
                [
                    atom.type + "_" + str(atom.residue.name)
                    for atom in self.structure_box_0_and_1_ff.atoms
                ]
            )

        else:
            self.types = np.array(
                [
                    atom.type + "_" + str(atom.residue.name)
                    for atom in self.structure_box_0_ff.atoms
                ]
            )

        self.unique_types = list(set(self.types))
        self.unique_types.sort(key=natural_sort)

        if self.structure_box_1:
            self.masses = (
                np.array(
                    [atom.mass for atom in self.structure_box_0_and_1_ff.atoms]
                )
                / self.mass_conversion_factor
            )
            self.mass_dict = dict(
                [
                    (self.unique_types.index(atom_type) + 1, mass)
                    for atom_type, mass in zip(self.types, self.masses)
                ]
            )

        else:
            self.masses = (
                np.array([atom.mass for atom in self.structure_box_0_ff.atoms])
                / self.mass_conversion_factor
            )
            self.mass_dict = dict(
                [
                    (self.unique_types.index(atom_type) + 1, mass)
                    for atom_type, mass in zip(self.types, self.masses)
                ]
            )

        # added an index so the atom types can be converted to numbers as the type name is to long for insertion into
        # the pdb and psf files
        self.atom_types_to_index_value_dict = dict(
            [
                (
                    self.unique_types[self.unique_types.index(atom_type)],
                    self.unique_types.index(atom_type),
                )
                for atom_type, mass in zip(self.types, self.masses)
            ]
        )

        # normalize by sigma
        self.box_0 = Box(
            lengths=np.array(
                [
                    (0.1 * val) / self.sigma_conversion_factor
                    for val in self.structure_box_0_ff.box[0:3]
                ]
            ),
            angles=self.structure_box_0_ff.box[3:6],
        )

        # create box 0 vector list and convert from nm to Ang and round to 6 decimals.
        # note mbuild standard lengths are in nm, so round to 6+1 = 7 then mutlipy by 10
        box_0_lengths_ang = (
            self.box_0.lengths[0] * 10,
            self.box_0.lengths[1] * 10,
            self.box_0.lengths[2] * 10,
        )
        self.box_0_vectors = _lengths_angles_to_vectors(
            box_0_lengths_ang, self.box_0.angles, precision=6
        )

        # Internally use nm
        if self.structure_box_1:
            self.box_1 = Box(
                lengths=np.array(
                    [
                        (0.1 * val) / self.sigma_conversion_factor
                        for val in self.structure_box_1_ff.box[0:3]
                    ]
                ),
                angles=self.structure_box_1_ff.box[3:6],
            )

            # create box 1 vector list and convert from nm to Ang and round to 6 decimals.
            # note mbuild standard lengths are in nm, so round to 6+1 = 7 then mutlipy by 10
            box_1_lengths_ang = (
                self.box_1.lengths[0] * 10,
                self.box_1.lengths[1] * 10,
                self.box_1.lengths[2] * 10,
            )
            self.box_1_vectors = _lengths_angles_to_vectors(
                box_1_lengths_ang, self.box_1.angles, precision=6
            )

        # if self.structure_box_1 != None:
        if self.structure_box_1:
            self.structure_selection = self.structure_box_0_and_1_ff
        else:
            self.structure_selection = self.structure_box_0_ff

        # Non-Bonded forces
        epsilons = (
            np.array([atom.epsilon for atom in self.structure_selection.atoms])
            / self.epsilon_conversion_factor
        )
        sigmas = (
            np.array([atom.sigma for atom in self.structure_selection.atoms])
            / self.sigma_conversion_factor
        )
        forcefields = [
            atom.type + "_" + atom.residue.name
            for atom in self.structure_selection.atoms
        ]
        residues_all_list = [
            atom.residue.name for atom in self.structure_selection.atoms
        ]

        self.epsilon_dict = dict(
            [
                (self.unique_types.index(atom_type), epsilon)
                for atom_type, epsilon in zip(self.types, epsilons)
            ]
        )
        self.sigma_dict = dict(
            [
                (self.unique_types.index(atom_type), sigma)
                for atom_type, sigma in zip(self.types, sigmas)
            ]
        )
        self.LJ_1_4_dict = dict(
            [
                (
                    self.unique_types.index(atom_type),
                    self.combined_1_4_lj_dict_per_residue[residues_all_list],
                )
                for atom_type, residues_all_list in zip(
                    self.types, residues_all_list
                )
            ]
        )
        self.forcefield_dict = dict(
            [
                (self.unique_types.index(atom_type), forcefield)
                for atom_type, forcefield in zip(self.types, forcefields)
            ]
        )

        # ensure all 1,4-coulombic scaling factors are the same
        coul_1_4_List = []
        for p in self.combined_1_4_coul_dict_per_residue.values():
            coul_1_4_List.append(p)
        coul_1_4_set = set(coul_1_4_List)
        if len(coul_1_4_set) > 1:
            self.input_error = True
            print_error_message = (
                "ERROR: There are multiple 1,4-coulombic scaling factors "
                "GOMC will only accept a singular input for the 1,4-coulombic "
                "scaling factors."
            )
            raise ValueError(print_error_message)
        else:
            self.coul_1_4 = list(coul_1_4_set)[0]

        # get all the unique atom name to check for the MEMC move in the gomc_conf_writer
        self.all_individual_atom_names_list = []
        self.all_residue_names_List = []
        if self.structure_box_1:
            list_of_structures = [
                self.structure_box_0_ff,
                self.structure_box_1_ff,
            ]
            stuct_only = [self.structure_box_0_ff, self.structure_box_1_ff]
        else:
            list_of_structures = [self.structure_box_0_ff]
            stuct_only = [self.structure_box_0_ff]

        for q_i in range(0, len(list_of_structures)):
            stuct_only_iteration = stuct_only[q_i]

            # caluculate the atom name and unique atom names
            residue_data_list = []
            residue_names_list = []
            for k, atom in enumerate(stuct_only_iteration.atoms):
                residue_data_list.append(str(atom.residue))
                residue_names_list.append(atom.residue.name)

            unique_residue_data_dict = {}
            unique_residue_data_list = []
            residue_data_name_list = []

            for m, residue in enumerate(stuct_only_iteration.residues):
                unique_residue_data_list.append(
                    str(stuct_only_iteration.residues[m])
                )
                unique_residue_data_dict.update(
                    {unique_residue_data_list[m]: m + 1}
                )
                residue_data_name_list.append(
                    stuct_only_iteration.residues[m].name
                )

            self.max_residue_no = 9999
            self.max_resname_char = 4

            res_no_chain_iter_corrected = []
            residue_id_list = []
            residue_id_adder_fixed_struct_wo_bonds = (
                0  # for example zeolite used as fixed atoms wo bonds
            )
            for f, PSF_atom_iteration_0 in enumerate(
                stuct_only_iteration.atoms
            ):

                if f > 0:
                    if (
                        PSF_atom_iteration_0.residue.chain
                        == previous_residue_chain
                        and len(PSF_atom_iteration_0.bonds) == 0
                    ):
                        residue_id_adder_fixed_struct_wo_bonds += 1

                previous_residue_chain = PSF_atom_iteration_0.residue.chain

                residue_id_int = int(
                    unique_residue_data_dict[residue_data_list[f]]
                    + residue_id_adder_fixed_struct_wo_bonds
                )

                res_id_adder = int(
                    (residue_id_int % self.max_residue_no) % self.max_residue_no
                )
                if int(res_id_adder) == 0:
                    res_no_iteration_corrected = int(self.max_residue_no)
                else:
                    res_no_iteration_corrected = res_id_adder

                res_no_chain_iter_corrected.append(res_no_iteration_corrected)
                residue_id_list.append(residue_id_int)

            # This converts the atom name in the GOMC psf and pdb files to unique atom names
            [
                unique_individual_atom_names_dict_iter,
                individual_atom_names_list_iter,
                missing_bead_to_atom_name_iter,
            ] = unique_atom_naming(
                stuct_only_iteration,
                residue_id_list,
                residue_names_list,
                bead_to_atom_name_dict=self.bead_to_atom_name_dict,
            )

            if q_i == 0:
                self.all_individual_atom_names_list = (
                    individual_atom_names_list_iter
                )

                self.all_residue_names_List = residue_names_list
            else:

                self.all_individual_atom_names_list = (
                    self.all_individual_atom_names_list
                    + individual_atom_names_list_iter
                )

                self.all_residue_names_List = (
                    self.all_residue_names_List + residue_names_list
                )

        # put the  self.all_individual_atom_names_list and self.all_residue_names_List in a list to match
        # the the atom name with a residue and find the unique matches
        if None in [
            unique_individual_atom_names_dict_iter,
            individual_atom_names_list_iter,
            missing_bead_to_atom_name_iter,
        ]:
            self.input_error = True
            print_error_message = (
                "ERROR: The unique_atom_naming function failed while "
                "running the charmm_writer function. Ensure the proper inputs are "
                "in the bead_to_atom_name_dict."
            )
            raise ValueError(print_error_message)

        else:
            self.all_res_unique_atom_name_dict = {}
            for res_i in range(0, len(self.all_individual_atom_names_list)):
                try:
                    if (
                        self.all_individual_atom_names_list[res_i]
                        not in self.all_res_unique_atom_name_dict[
                            self.all_residue_names_List[res_i]
                        ]
                    ):
                        self.all_res_unique_atom_name_dict.setdefault(
                            self.all_residue_names_List[res_i], []
                        ).append(self.all_individual_atom_names_list[res_i])
                except:
                    self.all_res_unique_atom_name_dict.setdefault(
                        self.all_residue_names_List[res_i], []
                    ).append(self.all_individual_atom_names_list[res_i])

        print(
            "all_res_unique_atom_name_dict = {}".format(
                self.all_res_unique_atom_name_dict
            )
        )

    def write_inp(self):
        """This write_inp function writes the Charmm style parameter (force field) file, which can be utilized
        in the GOMC and NAMD engines."""
        print("******************************")
        print("")
        print(
            "The charmm force field file writer (the write_inp function) is running"
        )

        if self.ff_filename is None:
            self.input_error = True
            print_error_message = (
                "ERROR: The force field file name was not specified and in the "
                "Charmm object. "
                "Therefore, the force field file (.inp) can not be written. "
                "Please use the force field file name when building the Charmm object, "
                "then use the write_inp function."
            )
            raise TypeError(print_error_message)
        else:

            print("******************************")
            print("")
            print(
                "The charmm force field file writer (the write_inp function) is running"
            )
            print("******************************")
            print("")
            print("writing the GOMC force field file ")
            date_time = datetime.datetime.today()

            unique_residue_data_dict = {}
            unique_residue_data_list = []
            residue_data_name_list = []

            # if self.structure_box_1 != None:
            if self.structure_box_1:
                residue_iterate = 0
                for m, residue in enumerate(self.structure_box_0_ff.residues):
                    residue_iterate = residue_iterate + 1
                    unique_residue_data_list.append(
                        str(self.structure_box_0_ff.residues[m])
                    )
                    unique_residue_data_dict.update(
                        {unique_residue_data_list[m]: m + 1}
                    )
                    residue_data_name_list.append(
                        self.structure_box_0_ff.residues[m].name
                    )

                for m, residue in enumerate(self.structure_box_1_ff.residues):
                    unique_residue_data_list.append(
                        str(self.structure_box_1_ff.residues[m])
                    )
                    unique_residue_data_dict.update(
                        {unique_residue_data_list[m]: m + 1 + residue_iterate}
                    )
                    residue_data_name_list.append(
                        self.structure_box_1_ff.residues[m].name
                    )

            else:
                for m, residue in enumerate(self.structure_box_0_ff.residues):
                    unique_residue_data_list.append(
                        str(self.structure_box_0_ff.residues[m])
                    )
                    unique_residue_data_dict.update(
                        {unique_residue_data_list[m]: m + 1}
                    )
                    residue_data_name_list.append(
                        self.structure_box_0_ff.residues[m].name
                    )

            for n in range(0, len(residue_data_name_list)):
                if residue_data_name_list[n] not in self.residues:
                    print(
                        "residue_data_name_list = "
                        + str(residue_data_name_list)
                    )
                    self.input_error = True
                    print_error_message = "ERROR: Please specifiy all residues (residues) in a list"
                    raise ValueError(print_error_message)

            with open(self.ff_filename, "w") as data:
                # Syntax which can change based on the functional form
                # Infer functional form based on the properties of the residue_id_list or structure_box_0_ff
                if self.detect_forcefield_style:
                    # Check angles
                    if len(self.structure_selection.urey_bradleys) > 0:
                        print("Urey bradley terms detected")
                        data.write(
                            "! Urey bradley terms detected, will use angle_style self"
                        )
                        data.write(
                            "! POTENTIAL ERROR: GOMC does no support the Urey bradley terms"
                        )
                        use_urey_bradleys = True
                    else:
                        print(
                            "No urey bradley terms detected, will use angle_style harmonic"
                        )
                        use_urey_bradleys = False

                    # Check dihedrals
                    if len(self.structure_selection.rb_torsions) > 0:
                        print(
                            "will use CHARMM_torsions  =  K0 + K1 * (1 + Cos[n1*(t) - (d1)] ) + "
                            + "K2 * (1 + Cos[n2*(t) - (d2)] ) + K3 * (1 + Cos[n3*(t) - (d3)] ) + "
                            + "K4 * (1 + Cos[n4*(t) - (d4)] ) + K5 * (1 + Cos[n5*(t) - (d5)] ) "
                        )
                        use_rb_torsions = True

                    else:
                        use_rb_torsions = False

                    if len(self.structure_selection.dihedrals) > 0:
                        print(
                            "Charmm dihedrals detected, will use dihedral_style charmm"
                        )
                        # this will need tested with a standard charmm input format before releasing it
                        use_dihedrals = True
                        self.input_error = True
                        print_error_message = (
                            "ERROR: use_dihedrals = {} "
                            "Charmm dihedrals not yet supported.".format(
                                use_dihedrals
                            )
                        )
                        raise ValueError(print_error_message)
                    else:
                        use_dihedrals = False

                if use_rb_torsions and use_dihedrals:
                    warn(
                        "Multiple dihedral styles detected, check your Forcefield XML and structure_selection"
                    )

                # Check impropers
                for dihedral in self.structure_selection.dihedrals:
                    if dihedral.improper:
                        warn(
                            "Amber-style impropers are currently not supported"
                        )

                bonds = [
                    [bond.atom1.idx + 1, bond.atom2.idx + 1]
                    for bond in self.structure_selection.bonds
                ]
                angles = [
                    [
                        angle.atom1.idx + 1,
                        angle.atom2.idx + 1,
                        angle.atom3.idx + 1,
                    ]
                    for angle in self.structure_selection.angles
                ]
                if use_rb_torsions:
                    dihedrals = [
                        [
                            dihedral.atom1.idx + 1,
                            dihedral.atom2.idx + 1,
                            dihedral.atom3.idx + 1,
                            dihedral.atom4.idx + 1,
                        ]
                        for dihedral in self.structure_selection.rb_torsions
                    ]
                elif use_dihedrals:
                    dihedrals = [
                        [
                            dihedral.atom1.idx + 1,
                            dihedral.atom2.idx + 1,
                            dihedral.atom3.idx + 1,
                            dihedral.atom4.idx + 1,
                        ]
                        for dihedral in self.structure_selection.dihedrals
                    ]
                else:
                    dihedrals = []
                impropers = [
                    [
                        improper.atom1.idx + 1,
                        improper.atom2.idx + 1,
                        improper.atom3.idx + 1,
                        improper.atom4.idx + 1,
                    ]
                    for improper in self.structure_selection.impropers
                ]

                if bonds:
                    if len(self.structure_selection.bond_types) == 0:
                        bond_types = np.ones(len(bonds), dtype=int)
                    else:
                        bond_types, unique_bond_types = _get_bond_types(
                            self.structure_selection,
                            self.sigma_conversion_factor,
                            self.epsilon_conversion_factor,
                        )

                if angles:
                    angle_types, unique_angle_types = _get_angle_types(
                        self.structure_selection,
                        self.sigma_conversion_factor,
                        self.epsilon_conversion_factor,
                        use_urey_bradleys=use_urey_bradleys,
                    )

                if dihedrals:
                    dihedral_types, unique_dihedral_types = _get_dihedral_types(
                        self.structure_selection,
                        use_rb_torsions,
                        use_dihedrals,
                        self.epsilon_conversion_factor,
                    )

                if impropers:
                    improper_types, unique_improper_types = _get_impropers(
                        self.structure_selection, self.epsilon_conversion_factor
                    )

                # if self.structure_box_1 != None:
                if self.structure_box_1:
                    data.write(
                        "*  "
                        + self.filename_box_0
                        + " and "
                        + self.filename_box_1
                        + " - created by mBuild using the on "
                        + str(date_time)
                        + "\n"
                    )
                else:
                    data.write(
                        "*  "
                        + self.filename_box_0
                        + " - created by mBuild using the on "
                        + str(date_time)
                        + "\n"
                    )

                data.write(
                    "*  "
                    + "parameters from the "
                    + str(self.forcefield_selection)
                    + " force field(s) via MoSDef\n"
                )
                data.write(
                    "*  1-4 coulombic scaling = "
                    + str(self.combined_1_4_coul_dict_per_residue)
                    + ", and 1-4 LJ scaling = "
                    + str(self.combined_1_4_lj_dict_per_residue)
                    + "\n\n"
                )
                data.write(
                    "*  "
                    + "{:d} atoms\n".format(len(self.structure_selection.atoms))
                )

                if self.atom_style in ["full", "molecular"]:
                    data.write("*  " + "{:d} bonds\n".format(len(bonds)))
                    data.write("*  " + "{:d} angles\n".format(len(angles)))
                    data.write(
                        "*  " + "{:d} dihedrals\n".format(len(dihedrals))
                    )
                    data.write(
                        "*  " + "{:d} impropers\n\n".format(len(impropers))
                    )

                data.write(
                    "*  " + "{:d} atom types\n".format(len(set(self.types)))
                )
                if self.atom_style in ["full", "molecular"]:
                    if bonds:
                        data.write(
                            "*  "
                            + "{:d} bond types\n".format(
                                len(set(unique_bond_types))
                            )
                        )
                    if angles:
                        data.write(
                            "*  "
                            + "{:d} angle types\n".format(
                                len(set(unique_angle_types))
                            )
                        )
                    if dihedrals:
                        data.write(
                            "*  "
                            + "{:d} dihedral types\n".format(
                                len(set(unique_dihedral_types))
                            )
                        )
                    if impropers:
                        data.write(
                            "*  "
                            + "{:d} improper types\n".format(
                                len(set(unique_improper_types))
                            )
                        )

                data.write("\n")

                data.write("\n* masses\n\n")
                data.write(
                    "!atom_types \tmass \t\t  atomTypeForceFieldName_ResidueName "
                    + "(i.e., atoms_type_per_utilized_FF)  \n"
                )
                for atom_type, mass in self.mass_dict.items():
                    mass_format = "*  {}\t\t{:.6f}\t! {}\n"
                    data.write(
                        mass_format.format(
                            base10_to_base52_alph(
                                atom_type - self.sub_1_for_base_52
                            ),
                            mass,
                            self.unique_types[atom_type - 1],
                        )
                    )

                # Bond coefficients
                if bonds:
                    data.write("\n")
                    data.write("BONDS * harmonic\n")
                    data.write("!\n")
                    data.write("!V(bond) = Kb(b - b0)**2\n")
                    data.write("!\n")
                    data.write("!Kb: kcal/mole/A**2\n")
                    data.write("!b0: A\n")
                    data.write(
                        "!Kb (kcal/mol) = Kb_K (K) * Boltz. const.; (9999999999 if no stretching)\n"
                    )
                    data.write("!\n")

                    if self.unit_style == "real":
                        data.write(
                            "!atom_types \t Kb\t\tb0 \t\t  atoms_types_per_utilized_FF\n"
                        )
                    for params, idx in unique_bond_types.items():
                        bond_format = "{}\t{}\t{}\t{}\t\t! {}\t{}\n"
                        if (
                            (self.gomc_fix_bonds_angles is not None)
                            and (
                                (params[3] and params[4])
                                in self.gomc_fix_bonds_angles
                            )
                        ) or (
                            (
                                (self.gomc_fix_bonds is not None)
                                and (
                                    (params[3] and params[4])
                                    in self.gomc_fix_bonds
                                )
                            )
                        ):
                            fix_bond_k_value = "999999999999"
                            data.write(
                                bond_format.format(
                                    base10_to_base52_alph(
                                        self.atom_types_to_index_value_dict[
                                            params[2][0] + "_" + str(params[3])
                                        ]
                                    ),
                                    base10_to_base52_alph(
                                        self.atom_types_to_index_value_dict[
                                            params[2][1] + "_" + str(params[4])
                                        ]
                                    ),
                                    fix_bond_k_value,
                                    params[1],
                                    params[2][0] + "_" + str(params[3]),
                                    params[2][1] + "_" + str(params[4]),
                                )
                            )

                        else:
                            data.write(
                                bond_format.format(
                                    base10_to_base52_alph(
                                        self.atom_types_to_index_value_dict[
                                            params[2][0] + "_" + str(params[3])
                                        ]
                                    ),
                                    base10_to_base52_alph(
                                        self.atom_types_to_index_value_dict[
                                            params[2][1] + "_" + str(params[4])
                                        ]
                                    ),
                                    params[0],
                                    params[1],
                                    params[2][0] + "_" + str(params[3]),
                                    params[2][1] + "_" + str(params[4]),
                                )
                            )

                # Angle coefficients
                if angles:
                    if use_urey_bradleys:
                        data.write(
                            "\n!  Urey Bradley terms detected but not written,"
                            + "since they are currently not compatible with GOMC\n"
                        )

                    data.write("\nANGLES * harmonic\n")
                    data.write("!\n")
                    data.write("!V(angle) = Ktheta(Theta - Theta0)**2\n")
                    data.write("!\n")
                    data.write("!Ktheta: kcal/mole/rad**2\n")
                    data.write("!Theta0: degrees\n")
                    data.write("!\n")
                    data.write(
                        "! Ktheta (kcal/mol) = Ktheta_K (K) * Boltz. const.\t\t\n"
                    )
                    data.write("!\n")
                    data.write(
                        "!atom_types \t\tKtheta\t\tTheta0\t\t\t  atoms_types_per_utilized_FF\n"
                    )
                    for params, idx in unique_angle_types.items():

                        if (
                            (self.gomc_fix_bonds_angles is not None)
                            and ((params[4] and params[5] and params[6]))
                            in self.gomc_fix_bonds_angles
                        ) or (
                            (self.gomc_fix_angles is not None)
                            and ((params[4] and params[5] and params[6]))
                            in self.gomc_fix_angles
                        ):
                            fix_angle_k_value = "999999999999"
                            angle_format = (
                                "{}\t{}\t{}\t{}\t\t{:.5f}\t\t! {}\t{}\t{}\n"
                            )
                            data.write(
                                angle_format.format(
                                    base10_to_base52_alph(
                                        self.atom_types_to_index_value_dict[
                                            params[3][0] + "_" + params[4]
                                        ]
                                    ),
                                    base10_to_base52_alph(
                                        self.atom_types_to_index_value_dict[
                                            params[2] + "_" + params[5]
                                        ]
                                    ),
                                    base10_to_base52_alph(
                                        self.atom_types_to_index_value_dict[
                                            params[3][1] + "_" + params[6]
                                        ]
                                    ),
                                    fix_angle_k_value,
                                    params[1],
                                    params[3][0] + "_" + params[4],
                                    params[2] + "_" + params[5],
                                    params[3][1] + "_" + params[6],
                                )
                            )

                        else:
                            data.write(
                                "{}\t{}\t{}\t{}\t\t{:.5f}\t\t! {}\t{}\t{}\n".format(
                                    base10_to_base52_alph(
                                        self.atom_types_to_index_value_dict[
                                            params[3][0] + "_" + params[4]
                                        ]
                                    ),
                                    base10_to_base52_alph(
                                        self.atom_types_to_index_value_dict[
                                            params[2] + "_" + params[5]
                                        ]
                                    ),
                                    base10_to_base52_alph(
                                        self.atom_types_to_index_value_dict[
                                            params[3][1] + "_" + params[6]
                                        ]
                                    ),
                                    params[0],
                                    params[1],
                                    params[3][0] + "_" + params[4],
                                    params[2] + "_" + params[5],
                                    params[3][1] + "_" + params[6],
                                )
                            )

                # Dihedral coefficients
                if dihedrals:
                    if use_rb_torsions:
                        list_if_large_error_dihedral_overall = []

                        list_if_largest_error_abs_values_for_dihedral_overall = (
                            []
                        )
                        list_dihedral_overall_error = []

                        list_if_abs_max_values_for_dihedral_overall = []
                        list_dihedral_atoms_all_dihedral_overall = []

                        data.write("\nDIHEDRALS * CHARMM\n")
                        data.write("!\n")
                        data.write(
                            "!V(dihedral) = Kchi(1 + cos(n(chi) - delta))\n"
                        )
                        data.write("!\n")
                        data.write("!Kchi: kcal/mole\n")
                        data.write("!n: multiplicity\n")
                        data.write("!delta: degrees\n")
                        data.write("!\n")
                        data.write(
                            "! Kchi (kcal/mol) = Kchi_K (K) * Boltz. const.\n"
                        )
                        data.write(
                            "! Boltzmann = 0.0019872041 kcal / (mol * K)\n"
                        )
                        data.write("!\n")
                        if self.unit_style == "real":
                            data.write(
                                "!atom_types \t\t\tKchi\t\tn\tdelta\t\t  atoms_types_per_utilized_FF\n"
                            )

                        for params, idx in unique_dihedral_types.items():
                            charmm_coeffs = RB_to_CHARMM(
                                params[0],
                                params[1],
                                params[2],
                                params[3],
                                params[4],
                                params[5],
                            )

                            # check the error between the convertions of RB_tortions to CHARMM DIHEDRALS (start)
                            rb_to_charmm_abs_diff = []
                            no_pi = np.pi
                            dihedral_steps = 2 * 10 ** (-3)
                            dihedral_range = 4 * no_pi
                            dihedral_no_steps = (
                                int(dihedral_range / dihedral_steps) + 1
                            )

                            for i in range(0, dihedral_no_steps + 1):
                                t = i * dihedral_steps

                                rb_dihedral_calc = (
                                    params[0]
                                    + params[1] * (np.cos(t - no_pi)) ** 1
                                    + params[2] * (np.cos(t - no_pi)) ** 2
                                    + params[3] * (np.cos(t - no_pi)) ** 3
                                    + params[4] * (np.cos(t - no_pi)) ** 4
                                    + params[5] * (np.cos(t - no_pi)) ** 5
                                )
                                """CHARMM_torsions
                                = K0 * (1 + Cos[n0 * (t) - (d0)]) + K1 * (1 + Cos[n1 * (t) - (d1)]) + K2 * (
                                           # 1 + Cos[n2 * (t) - (d2)])
                                + K3 * (1 + Cos[n3 * (t) - (d3)]) + K4 * (1 + Cos[n4 * (t) - (d4)]) + K5 * (
                                           # 1 + Cos[n5 * (t) - (d5)])

                                = K0 + K1 * (1 + Cos[n1 * (t) - (d1)]) + K2 * (1 + Cos[n2 * (t) - (d2)])
                                + K3 * (1 + Cos[n3 * (t) - (d3)]) + K4 * (1 + Cos[n4 * (t) - (d4)]) + K5 * (
                                           # 1 + Cos[n5 * (t) - (d5)]). """

                                rb_to_charmm_calc = (
                                    charmm_coeffs[0, 0]
                                    * (
                                        1
                                        + np.cos(
                                            int(charmm_coeffs[0, 1]) * t
                                            - charmm_coeffs[0, 2] * no_pi / 180
                                        )
                                    )
                                    + charmm_coeffs[1, 0]
                                    * (
                                        1
                                        + np.cos(
                                            int(charmm_coeffs[1, 1]) * t
                                            - charmm_coeffs[1, 2] * no_pi / 180
                                        )
                                    )
                                    + charmm_coeffs[2, 0]
                                    * (
                                        1
                                        + np.cos(
                                            int(charmm_coeffs[2, 1]) * t
                                            - charmm_coeffs[2, 2] * no_pi / 180
                                        )
                                    )
                                    + charmm_coeffs[3, 0]
                                    * (
                                        1
                                        + np.cos(
                                            int(charmm_coeffs[3, 1]) * t
                                            - charmm_coeffs[3, 2] * no_pi / 180
                                        )
                                    )
                                    + charmm_coeffs[4, 0]
                                    * (
                                        1
                                        + np.cos(
                                            int(charmm_coeffs[4, 1]) * t
                                            - charmm_coeffs[4, 2] * no_pi / 180
                                        )
                                    )
                                    + charmm_coeffs[5, 0]
                                    * (
                                        1
                                        + np.cos(
                                            int(charmm_coeffs[5, 1]) * t
                                            - charmm_coeffs[5, 2] * no_pi / 180
                                        )
                                    )
                                )

                                rb_to_charmm_absolute_difference = np.absolute(
                                    rb_dihedral_calc - rb_to_charmm_calc
                                )
                                rb_to_charmm_abs_diff.append(
                                    rb_to_charmm_absolute_difference
                                )

                            list_if_large_error_dihedral_iteration = []
                            list_abs_max_dihedral_iteration = []

                            if max(rb_to_charmm_abs_diff) > 10 ** (-10):
                                list_if_large_error_dihedral_iteration.append(1)
                                list_abs_max_dihedral_iteration.append(
                                    max(rb_to_charmm_abs_diff)
                                )

                                list_if_large_error_dihedral_overall.append(1)
                                list_if_largest_error_abs_values_for_dihedral_overall.append(
                                    max(rb_to_charmm_abs_diff)
                                )
                                list_dihedral_overall_error.append(
                                    str(params[8])
                                    + ", "
                                    + str(params[9])
                                    + ", "
                                    + str(params[10])
                                    + ", "
                                    + str(params[11])
                                )

                            else:
                                list_if_large_error_dihedral_iteration.append(0)

                                list_if_abs_max_values_for_dihedral_overall.append(
                                    max(rb_to_charmm_abs_diff)
                                )
                                list_dihedral_atoms_all_dihedral_overall.append(
                                    str(params[8])
                                    + ", "
                                    + str(params[9])
                                    + ", "
                                    + str(params[10])
                                    + ", "
                                    + str(params[11])
                                )

                            # **************************************
                            # check the error between the convertions of RB_tortions to CHARMM DIHEDRALS (end)
                            # **************************************
                            dihedral_format = "{}\t{}\t{}\t{}\t{:.6f}\t{}\t{}\t\t! {}\t{}\t{}\t{}\n"

                            # Note the Charmm C0 or K0 dihedral term is not printed as it is used/read
                            # as a harmonic potential in the Charmm format

                            data.write(
                                dihedral_format.format(
                                    base10_to_base52_alph(
                                        self.atom_types_to_index_value_dict[
                                            params[8] + "_" + params[12]
                                        ]
                                    ),
                                    base10_to_base52_alph(
                                        self.atom_types_to_index_value_dict[
                                            params[9] + "_" + params[13]
                                        ]
                                    ),
                                    base10_to_base52_alph(
                                        self.atom_types_to_index_value_dict[
                                            params[10] + "_" + params[14]
                                        ]
                                    ),
                                    base10_to_base52_alph(
                                        self.atom_types_to_index_value_dict[
                                            params[11] + "_" + params[15]
                                        ]
                                    ),
                                    charmm_coeffs[1, 0],
                                    int(charmm_coeffs[1, 1]),
                                    charmm_coeffs[1, 2],
                                    params[8] + "_" + params[12],
                                    params[9] + "_" + params[13],
                                    params[10] + "_" + params[14],
                                    params[11] + "_" + params[15],
                                )
                            )
                            data.write(
                                dihedral_format.format(
                                    base10_to_base52_alph(
                                        self.atom_types_to_index_value_dict[
                                            params[8] + "_" + params[12]
                                        ]
                                    ),
                                    base10_to_base52_alph(
                                        self.atom_types_to_index_value_dict[
                                            params[9] + "_" + params[13]
                                        ]
                                    ),
                                    base10_to_base52_alph(
                                        self.atom_types_to_index_value_dict[
                                            params[10] + "_" + params[14]
                                        ]
                                    ),
                                    base10_to_base52_alph(
                                        self.atom_types_to_index_value_dict[
                                            params[11] + "_" + params[15]
                                        ]
                                    ),
                                    charmm_coeffs[2, 0],
                                    int(charmm_coeffs[2, 1]),
                                    charmm_coeffs[2, 2],
                                    params[8] + "_" + params[12],
                                    params[9] + "_" + params[13],
                                    params[10] + "_" + params[14],
                                    params[11] + "_" + params[15],
                                )
                            )
                            data.write(
                                dihedral_format.format(
                                    base10_to_base52_alph(
                                        self.atom_types_to_index_value_dict[
                                            params[8] + "_" + params[12]
                                        ]
                                    ),
                                    base10_to_base52_alph(
                                        self.atom_types_to_index_value_dict[
                                            params[9] + "_" + params[13]
                                        ]
                                    ),
                                    base10_to_base52_alph(
                                        self.atom_types_to_index_value_dict[
                                            params[10] + "_" + params[14]
                                        ]
                                    ),
                                    base10_to_base52_alph(
                                        self.atom_types_to_index_value_dict[
                                            params[11] + "_" + params[15]
                                        ]
                                    ),
                                    charmm_coeffs[3, 0],
                                    int(charmm_coeffs[3, 1]),
                                    charmm_coeffs[3, 2],
                                    params[8] + "_" + params[12],
                                    params[9] + "_" + params[13],
                                    params[10] + "_" + params[14],
                                    params[11] + "_" + params[15],
                                )
                            )
                            data.write(
                                dihedral_format.format(
                                    base10_to_base52_alph(
                                        self.atom_types_to_index_value_dict[
                                            params[8] + "_" + params[12]
                                        ]
                                    ),
                                    base10_to_base52_alph(
                                        self.atom_types_to_index_value_dict[
                                            params[9] + "_" + params[13]
                                        ]
                                    ),
                                    base10_to_base52_alph(
                                        self.atom_types_to_index_value_dict[
                                            params[10] + "_" + params[14]
                                        ]
                                    ),
                                    base10_to_base52_alph(
                                        self.atom_types_to_index_value_dict[
                                            params[11] + "_" + params[15]
                                        ]
                                    ),
                                    charmm_coeffs[4, 0],
                                    int(charmm_coeffs[4, 1]),
                                    charmm_coeffs[4, 2],
                                    params[8] + "_" + params[12],
                                    params[9] + "_" + params[13],
                                    params[10] + "_" + params[14],
                                    params[11] + "_" + params[15],
                                )
                            )
                            data.write(
                                dihedral_format.format(
                                    base10_to_base52_alph(
                                        self.atom_types_to_index_value_dict[
                                            params[8] + "_" + params[12]
                                        ]
                                    ),
                                    base10_to_base52_alph(
                                        self.atom_types_to_index_value_dict[
                                            params[9] + "_" + params[13]
                                        ]
                                    ),
                                    base10_to_base52_alph(
                                        self.atom_types_to_index_value_dict[
                                            params[10] + "_" + params[14]
                                        ]
                                    ),
                                    base10_to_base52_alph(
                                        self.atom_types_to_index_value_dict[
                                            params[11] + "_" + params[15]
                                        ]
                                    ),
                                    charmm_coeffs[5, 0],
                                    int(charmm_coeffs[5, 1]),
                                    charmm_coeffs[5, 2],
                                    params[8] + "_" + params[12],
                                    params[9] + "_" + params[13],
                                    params[10] + "_" + params[14],
                                    params[11] + "_" + params[15],
                                )
                            )

                        if sum(list_if_large_error_dihedral_overall) > 0:
                            list_max_error_abs_dihedral_overall = max(
                                list_if_largest_error_abs_values_for_dihedral_overall
                            )
                            info_if_dihedral_error_too_large = (
                                "! WARNING: RB-torsion to CHARMM "
                                "dihedral conversion error"
                                " is to large [error > 10^(-10)] \n"
                                "! WARNING: Maximum( "
                                "|(RB-torsion calc)-(CHARMM dihedral calc)| ) =  "
                                + str(list_max_error_abs_dihedral_overall)
                                + "\n"
                            )
                            warn(
                                "! WARNING: RB-torsion to CHARMM dihedral conversion error"
                                "is to large [error > 10^(-10)] \n"
                                "! WARNING: Maximum( |(RB-torsion calc)-(CHARMM dihedral calc)| ) =  "
                                + str(list_max_error_abs_dihedral_overall)
                                + "\n"
                            )
                            data.write(info_if_dihedral_error_too_large)
                            print(info_if_dihedral_error_too_large)
                        else:
                            list_if_abs_max_values_for_dihedral_overall_max = (
                                max(list_if_abs_max_values_for_dihedral_overall)
                            )
                            info_if_dihedral_error_ok = (
                                "! RB-torsion to CHARMM dihedral conversion error is OK "
                                "[error <= 10^(-10)]\n"
                                + "! Maximum( |(RB-torsion calc)-(CHARMM dihedral calc)| ) =  "
                                + str(
                                    list_if_abs_max_values_for_dihedral_overall_max
                                )
                                + "\n"
                            )
                            data.write(info_if_dihedral_error_ok)
                            print(info_if_dihedral_error_ok)

                    elif use_dihedrals:
                        data.write(
                            "ERROR: not set up to use to use_dihedrals form for data input from the xml file"
                        )

                # Improper coefficients
                if impropers:
                    data.write(
                        "ERROR: GOMC is not currently able to use improper in its calculations"
                    )

                # Pair coefficients
                print(
                    "NBFIX_Mixing not used or no mixing used for the non-bonded potentials out"
                )
                print("self.non_bonded_type = " + str(self.non_bonded_type))
                if self.non_bonded_type == "LJ":
                    data.write("\n")
                    data.write("NONBONDED\n")
                    data.write("!\n")
                    data.write(
                        "!V(Lennard-Jones) = Eps,i,j[(Rmin,i,j/ri,j)**12 - 2(Rmin,i,j/ri,j)**6]\n"
                    )
                    data.write("!\n")

                    if self.unit_style == "real":
                        data.write(
                            "!atype \tignored\tepsilon \tRmin/2 \t\tignored\teps,1-4\t\tRmin/2,1-4\t\t"
                            + "  atom_type_per_utilized_FF\n"
                        )

                    print("forcefield_dict = " + str(self.forcefield_dict))

                    for idx, epsilon in self.epsilon_dict.items():
                        nb_format = "{}\t{:.2f}\t{:.9f}\t{:.11f}\t{:.2f}\t{:.9f}\t{:.11f}\t\t! {}\t{}\n"
                        data.write(
                            nb_format.format(
                                base10_to_base52_alph(idx),
                                0,
                                -epsilon,
                                self.sigma_dict[idx] * (2 ** (1 / 6)) / 2,
                                0,
                                float(self.LJ_1_4_dict[idx]) * (-epsilon),
                                self.sigma_dict[idx] * (2 ** (1 / 6)) / 2,
                                self.forcefield_dict[idx],
                                self.forcefield_dict[idx],
                            )
                        )

                elif self.non_bonded_type == "Mie":
                    data.write(
                        "ERROR: Currently the Mie potential (non_bonded_type) is not supported in this MoSDeF "
                        "GOMC parameter writer\n"
                    )
                    print_error_message = (
                        "ERROR: Currently the Mie potential (non_bonded_type) is not "
                        "supported in this MoSDeF GOMC parameter writer."
                    )
                    raise ValueError(print_error_message)
                else:
                    data.write(
                        "ERROR: Currently this potential (non_bonded_type) is not supported in "
                        "this MoSDeF GOMC parameter writer\n"
                    )
                    print_error_message = (
                        "ERROR: Currently this potential (non_bonded_type) is not supported in "
                        "this MoSDeF GOMC parameter writer."
                    )
                    raise ValueError(print_error_message)

                # writing end in file
                data.write("\nEND\n")

        # **********************************
        # **********************************
        # FF writer (end)
        # **********************************
        # **********************************

    def write_psf(self):
        """This write_psf function writes the Charmm style PSF (topology) file, which can be utilized
        in the GOMC and NAMD engines."""
        # **********************************
        # **********************************
        # psf writer (start)
        # **********************************
        # **********************************

        print("******************************")
        print("")
        print(
            "The charmm X-plor format psf writer (the write_psf function) is running"
        )

        date_time = datetime.datetime.today()

        print(
            "write_psf: forcefield_selection = {}, residues = {}".format(
                self.forcefield_selection, self.residues
            )
        )

        print("******************************")
        print("")

        if self.structure_box_1:
            list_of_structures = [
                self.structure_box_0_ff,
                self.structure_box_1_ff,
            ]
            list_of_file_names = [self.filename_box_0, self.filename_box_1]
            stuct_only = [self.structure_box_0_ff, self.structure_box_1_ff]
        else:
            list_of_structures = [self.structure_box_0_ff]
            list_of_file_names = [self.filename_box_0]
            stuct_only = [self.structure_box_0_ff]

        for q in range(0, len(list_of_structures)):
            stuct_iteration = list_of_structures[q]
            file_name_iteration = list_of_file_names[q]
            output = str(file_name_iteration) + ".psf"
            stuct_only_iteration = stuct_only[q]
            # Lammps syntax depends on the functional form
            # Infer functional form based on the properties of the stuct_iteration
            if self.detect_forcefield_style:
                # Check  for angles
                if len(stuct_iteration.urey_bradleys) > 0:
                    print(
                        "Warning: Urey bradley terms detected. GOMC does no support the Urey-Bradley terms"
                    )
                    warn(
                        "warning: Urey bradley terms detected. "
                        "GOMC does no support the Urey-Bradley terms"
                    )
                    use_urey_bradleys = True
                else:
                    print("No urey bradley terms detected")
                    use_urey_bradleys = False

                # Check for dihedrals
                if len(stuct_iteration.rb_torsions) > 0:
                    print(
                        "RB Torsions detected, will converted to CHARMM Dihedrals"
                    )
                    use_rb_torsions = True
                    dihedrals_list = stuct_iteration.rb_torsions
                    dihedrals = [
                        [
                            dihedral.atom1.idx + 1,
                            dihedral.atom2.idx + 1,
                            dihedral.atom3.idx + 1,
                            dihedral.atom4.idx + 1,
                        ]
                        for dihedral in stuct_iteration.rb_torsions
                    ]
                else:
                    use_rb_torsions = False

                if len(stuct_iteration.dihedrals) > 0:
                    print(
                        "Charmm dihedrals detected, so CHARMM Dihedrals will remain"
                    )
                    use_dihedrals = True
                    dihedrals_list = stuct_iteration.dihedrals
                    dihedrals = [
                        [
                            dihedral.atom1.idx + 1,
                            dihedral.atom2.idx + 1,
                            dihedral.atom3.idx + 1,
                            dihedral.atom4.idx + 1,
                        ]
                        for dihedral in stuct_iteration.dihedrals
                    ]
                else:
                    use_dihedrals = False
            if (use_rb_torsions is False) and (use_dihedrals is False):
                dihedrals_list = []
                dihedrals = []
            if use_rb_torsions and use_dihedrals:
                warn(
                    "Multiple dihedral styles detected, check your "
                    "Forcefield XML and structure files"
                )

            # Check for impropers
            for dihedral in stuct_iteration.dihedrals:
                if dihedral.improper:
                    warn(
                        "ERROR: Amber-style impropers are currently not supported in GOMC"
                    )

            impropers_list = stuct_iteration.impropers
            impropers = [
                [
                    improper.atom1.idx + 1,
                    improper.atom2.idx + 1,
                    improper.atom3.idx + 1,
                    improper.atom4.idx + 1,
                ]
                for improper in stuct_iteration.impropers
            ]

            no_atoms = len(stuct_iteration.atoms)
            no_bonds = len(stuct_iteration.bonds)
            no_angles = len(stuct_iteration.angles)

            no_dihedrals = len(dihedrals)
            no_impropers = len(impropers)

            no_donors = len(stuct_iteration.donors)
            no_acceptors = len(stuct_iteration.acceptors)
            no_groups = len(stuct_iteration.groups)

            # psf printing (start)

            residue_data_list = []
            residue_names_list = []
            for k, atom in enumerate(stuct_only_iteration.atoms):
                residue_data_list.append(str(atom.residue))
                residue_names_list.append(atom.residue.name)

            unique_residue_data_dict = {}
            unique_residue_data_list = []
            residue_data_name_list = []

            for m, residue in enumerate(stuct_only_iteration.residues):
                unique_residue_data_list.append(
                    str(stuct_only_iteration.residues[m])
                )
                unique_residue_data_dict.update(
                    {unique_residue_data_list[m]: m + 1}
                )
                residue_data_name_list.append(
                    stuct_only_iteration.residues[m].name
                )

            res_no_chain_iter_corrected = []
            residue_id_list = []
            residue_id_adder_fixed_struct_wo_bonds = 0
            for f, PSF_atom_iteration_0 in enumerate(
                stuct_only_iteration.atoms
            ):
                if f > 0:
                    if (
                        PSF_atom_iteration_0.residue.chain
                        == previous_residue_chain
                        and len(PSF_atom_iteration_0.bonds) == 0
                    ):
                        residue_id_adder_fixed_struct_wo_bonds += 1

                previous_residue_chain = PSF_atom_iteration_0.residue.chain

                residue_id_int = int(
                    unique_residue_data_dict[residue_data_list[f]]
                    + residue_id_adder_fixed_struct_wo_bonds
                )
                res_id_adder = int(
                    (residue_id_int % self.max_residue_no) % self.max_residue_no
                )
                if int(res_id_adder) == 0:
                    res_no_iteration_corrected = int(self.max_residue_no)
                else:
                    res_no_iteration_corrected = res_id_adder

                res_no_chain_iter_corrected.append(res_no_iteration_corrected)
                residue_id_list.append(residue_id_int)

            output_write = genopen(output, "w")

            first_indent = "%8s"
            psf_formating = (
                "%8s %-4s %-4s %-4s %-4s %4s %10.6f %13.4f" + 11 * " "
            )

            output_write.write("PSF ")
            output_write.write("\n\n")

            no_of_remarks = 3
            output_write.write(first_indent % no_of_remarks + " !NTITLE\n")
            output_write.write(
                " REMARKS this file "
                + file_name_iteration
                + " - created by mBuild/foyer using the"
                + "\n"
            )
            output_write.write(
                " REMARKS parameters from the "
                + str(self.forcefield_selection)
                + " force field via MoSDef\n"
            )
            output_write.write(
                " REMARKS created on " + str(date_time) + "\n\n\n"
            )

            # This converts the atom name in the GOMC psf and pdb files to unique atom names
            print(
                "bead_to_atom_name_dict = {}".format(
                    self.bead_to_atom_name_dict
                )
            )
            [
                unique_individual_atom_names_dict,
                individual_atom_names_list,
                missing_bead_to_atom_name,
            ] = unique_atom_naming(
                stuct_only_iteration,
                residue_id_list,
                residue_names_list,
                bead_to_atom_name_dict=self.bead_to_atom_name_dict,
            )

            if None in [
                unique_individual_atom_names_dict,
                individual_atom_names_list,
                missing_bead_to_atom_name,
            ]:
                self.input_error = True
                print_error_message = (
                    "ERROR: The unique_atom_naming function failed while "
                    "running the charmm_writer function. Ensure the proper inputs are "
                    "in the bead_to_atom_name_dict."
                )
                raise ValueError(print_error_message)

            # ATOMS: Calculate the atom data
            # psf_formating is conducted for the for CHARMM format (i.e., atom types are base 52, letters only)
            output_write.write(first_indent % no_atoms + " !NATOM\n")
            for i_atom, PSF_atom_iteration_1 in enumerate(
                stuct_iteration.atoms
            ):
                segment_id = PSF_atom_iteration_1.residue.segid or "SYS"
                atom_type_iter = base10_to_base52_alph(
                    self.atom_types_to_index_value_dict[
                        PSF_atom_iteration_1.type
                        + "_"
                        + PSF_atom_iteration_1.residue.name
                    ]
                )

                atom_lines_iteration = psf_formating % (
                    i_atom + 1,
                    segment_id,
                    res_no_chain_iter_corrected[i_atom],
                    str(residue_names_list[i_atom])[: self.max_resname_char],
                    individual_atom_names_list[i_atom],
                    atom_type_iter,
                    PSF_atom_iteration_1.charge,
                    PSF_atom_iteration_1.mass,
                )

                output_write.write("%s\n" % atom_lines_iteration)

            output_write.write("\n")

            # BONDS: Calculate the bonding data
            output_write.write(first_indent % no_bonds + " !NBOND: bonds\n")
            for i_bond, PSF_bond_iteration_1 in enumerate(
                stuct_iteration.bonds
            ):
                output_write.write(
                    (first_indent * 2)
                    % (
                        PSF_bond_iteration_1.atom1.idx + 1,
                        PSF_bond_iteration_1.atom2.idx + 1,
                    )
                )

                if (i_bond + 1) % 4 == 0:
                    output_write.write("\n")

            if no_bonds % 4 == 0:
                output_write.write("\n")
            else:
                output_write.write("\n\n")

            if no_bonds == 0:
                output_write.write("\n")

            # ANGLES: Calculate the angle data
            output_write.write(first_indent % no_angles + " !NTHETA: angles\n")
            for i_angle, angle_iteration in enumerate(stuct_iteration.angles):

                output_write.write(
                    (first_indent * 3)
                    % (
                        angle_iteration.atom1.idx + 1,
                        angle_iteration.atom2.idx + 1,
                        angle_iteration.atom3.idx + 1,
                    )
                )

                if (i_angle + 1) % 3 == 0:
                    output_write.write("\n")

            if no_angles % 3 == 0:
                output_write.write("\n")
            else:
                output_write.write("\n\n")

            if no_angles == 0:
                output_write.write("\n")

            # DIHEDRALS: Calculate the dihedral  data
            output_write.write(
                first_indent % no_dihedrals + " !NPHI: dihedrals\n"
            )
            for i_dihedral, dihedral_iter in enumerate(dihedrals_list):
                (
                    dihedral_atom_1,
                    dihedral_atom_2,
                    dihedral_atom_3,
                    dihedral_atom_4,
                ) = (
                    dihedral_iter.atom1,
                    dihedral_iter.atom2,
                    dihedral_iter.atom3,
                    dihedral_iter.atom4,
                )

                output_write.write(
                    (first_indent * 4)
                    % (
                        dihedral_atom_1.idx + 1,
                        dihedral_atom_2.idx + 1,
                        dihedral_atom_3.idx + 1,
                        dihedral_atom_4.idx + 1,
                    )
                )

                if (i_dihedral + 1) % 2 == 0:
                    output_write.write("\n")

            if no_dihedrals % 2 == 0:
                output_write.write("\n")
            else:
                output_write.write("\n\n")

            if no_dihedrals == 0:
                output_write.write("\n")

            # IMPROPERS: Calculate the improper data
            output_write.write(
                first_indent % no_impropers + " !NIMPHI: impropers\n"
            )
            for i_improper, improper_iter in enumerate(impropers_list):
                (
                    improper_atom_1,
                    improper_atom_2,
                    improper_atom_3,
                    improper_atom_4,
                ) = (
                    improper_iter.atom1,
                    improper_iter.atom2,
                    improper_iter.atom3,
                    improper_iter.atom4,
                )

                output_write.write(
                    (first_indent * 4)
                    % (
                        improper_atom_1.idx + 1,
                        improper_atom_2.idx + 1,
                        improper_atom_3.idx + 1,
                        improper_atom_4.idx + 1,
                    )
                )

                if (i_improper + 1) % 2 == 0:
                    output_write.write("\n")

            if no_impropers % 2 == 0:
                output_write.write("\n")
            else:
                output_write.write("\n\n")

            if no_impropers == 0:
                output_write.write("\n")

            # DONOR: calculate the donor data
            output_write.write(first_indent % no_donors + " !NDON: donors\n")
            for donor_i, donor_iter in enumerate(stuct_iteration.donors):

                output_write.write(
                    (first_indent * 2)
                    % (donor_iter.atom1.idx + 1, donor_iter.atom2.idx + 1)
                )
                if (donor_i + 1) % 4 == 0:
                    output_write.write("\n")

            if no_donors % 4 == 0:
                output_write.write("\n")
            else:
                output_write.write("\n\n")

            if no_donors == 0:
                output_write.write("\n")

            # ACCEPTOR: calculate the acceptor data
            output_write.write(
                first_indent % no_acceptors + " !NACC: acceptors\n"
            )
            for acceptor_i, acceptor_iter in enumerate(
                stuct_iteration.acceptors
            ):

                output_write.write(
                    (first_indent * 2)
                    % (acceptor_iter.atom1.idx + 1, acceptor_iter.atom2.idx + 1)
                )
                if (acceptor_i + 1) % 4 == 0:
                    output_write.write("\n")

            if no_acceptors % 4 == 0:
                output_write.write("\n")
            else:
                output_write.write("\n\n")

            if no_acceptors == 0:
                output_write.write("\n")

            # NNB: calculate the NNB data
            output_write.write(first_indent % 0 + " !NNB\n\n")
            for nbb_i, atoms_iter in enumerate(stuct_iteration.atoms):

                output_write.write(first_indent % 0)
                if (nbb_i + 1) % 8 == 0:
                    output_write.write("\n")

            if no_atoms % 8 == 0:
                output_write.write("\n")
            else:
                output_write.write("\n\n")

            if no_atoms == 0:
                output_write.write("\n")

            # GROUP: calculate the group data
            try:
                group_data = stuct_iteration.groups.nst2
            except AttributeError:
                group_data = 0
            output_write.write(
                (first_indent * 2) % (no_groups or 1, group_data) + " !NGRP \n"
            )
            if stuct_iteration.groups is True:
                for group_i, group_iter in enumerate(stuct_iteration.groups):

                    output_write.write(
                        (first_indent * 3)
                        % (
                            group_iter.atom.idx,
                            group_iter.type,
                            group_iter.move,
                        )
                    )
                    if (group_i + 1) % 3 == 0:
                        output_write.write("\n")

                if no_groups % 3 == 0:
                    output_write.write("\n")
                else:
                    output_write.write("\n\n")

                if no_groups == 0:
                    output_write.write("\n")

            else:
                structure_abs_charge_value = abs(
                    sum(
                        atom_charge_iter.charge
                        for atom_charge_iter in stuct_iteration.atoms
                    )
                )
                if structure_abs_charge_value < 1.0e-4:
                    group_type = 1
                else:
                    group_type = 2
                output_write.write((first_indent * 3) % (0, group_type, 0))
                output_write.write("\n")

            output_write.write("\n")
            output_write.close()
        # **********************************
        # **********************************
        # psf writer (end)
        # **********************************
        # **********************************

    def write_pdb(self):
        """This write_psf function writes the Charmm style PDB (coordinate file), which can be utilized
        in the GOMC and NAMD engines."""
        # **********************************
        # **********************************
        # pdb writer (start)
        # **********************************
        # **********************************
        date_time = datetime.datetime.today()
        print("******************************")
        print("")
        print("The charmm pdb writer (the write_pdb function) is running")
        print("write_charmm_pdb: residues == {}".format(self.residues))
        print("fix_residue = {}".format(self.fix_residue))
        print("fix_residue_in_box = {}".format(self.fix_residue_in_box))
        print("bead_to_atom_name_dict = {}".format(self.bead_to_atom_name_dict))

        if self.fix_residue is None and self.fix_residue_in_box is None:
            print(
                "INFORMATION: No atoms are fixed in this pdb file for the GOMC simulation engine. "
            )
        else:
            warn(
                "Some atoms are fixed in this pdb file for the GOMC simulation engine. "
            )

        print("******************************")
        print("")

        if self.structure_box_1:
            list_of_structures = [
                self.structure_box_0_ff,
                self.structure_box_1_ff,
            ]
            list_of_file_names = [self.filename_box_0, self.filename_box_1]
            stuct_only = [self.structure_box_0_ff, self.structure_box_1_ff]
        else:
            list_of_structures = [self.structure_box_0_ff]
            list_of_file_names = [self.filename_box_0]
            stuct_only = [self.structure_box_0_ff]

        for q in range(0, len(list_of_structures)):
            file_name_iteration = list_of_file_names[q]
            output = str(file_name_iteration) + ".pdb"
            stuct_only_iteration = stuct_only[q]

            output_write = genopen(output, "w")
            # output_write.write(
            #'REMARK this file ' + file_name_iteration + ' - created by mBuild/foyer using the' + '\n')
            # output_write.write(
            #'REMARK parameters from the ' + str(self.forcefield_selection) + ' force field via MoSDef\n')
            # output_write.write('REMARK created on ' + str(date_time) + '\n')

            unique_residue_data_dict = {}
            unique_residue_data_list = []
            residue_data_name_list = []
            for m, residue in enumerate(stuct_only_iteration.residues):
                unique_residue_data_list.append(
                    str(stuct_only_iteration.residues[m])
                )
                unique_residue_data_dict.update(
                    {unique_residue_data_list[m]: m + 1}
                )
                residue_data_name_list.append(
                    stuct_only_iteration.residues[m].name
                )

            for n in range(0, len(residue_data_name_list)):
                if residue_data_name_list[n] not in self.residues:
                    self.input_error = True
                    print_error_message = "ERROR: Please specifiy all residues (residues) in a list"
                    raise ValueError(print_error_message)

            residue_data_list = []
            for k, atom in enumerate(stuct_only_iteration.atoms):
                residue_data_list.append(str(atom.residue))

            if (self.fix_residue is not None) and (
                self.fix_residue_in_box is not None
            ):
                for n in range(0, len(self.fix_residue)):
                    if self.fix_residue[n] in self.fix_residue_in_box:
                        self.input_error = True
                        print_error_message = (
                            "ERROR: residue type can not be specified to both "
                            "fix_residue and fix_residue_in_box"
                        )
                        raise ValueError(print_error_message)

            residue_names_list = []
            fix_atoms_list = []
            for k, atom in enumerate(stuct_only_iteration.atoms):
                residue_names_list.append(atom.residue.name)
                if (self.fix_residue is not None) and (
                    atom.residue.name in self.fix_residue
                ):
                    beta_iteration = 1.00
                elif (self.fix_residue_in_box is not None) and (
                    atom.residue.name in self.fix_residue_in_box
                ):
                    beta_iteration = 2.00
                else:
                    beta_iteration = 0.00
                fix_atoms_list.append(beta_iteration)

            if stuct_only_iteration.box is not None:
                output_write.write(
                    "CRYST1%9.3f%9.3f%9.3f%7.2f%7."
                    "2f%7.2f %-11s%4s\n"
                    % (
                        stuct_only_iteration.box[0],
                        stuct_only_iteration.box[1],
                        stuct_only_iteration.box[2],
                        stuct_only_iteration.box[3],
                        stuct_only_iteration.box[4],
                        stuct_only_iteration.box[5],
                        stuct_only_iteration.space_group,
                        "",
                    )
                )

            all_atom_coordinates = stuct_only_iteration.get_coordinates("all")
            if all_atom_coordinates is None:
                self.input_error = True
                print_error_message = (
                    "ERROR: the submitted structure has no PDB coordinates, "
                    "so the PDB writer has terminated. "
                )
                raise ValueError(print_error_message)

            pdb_atom_line_format = "ATOM  %5s %-4s%1s%-4s%1s%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f      %-4s%2s%-2s\n"

            atom_alternate_location_list = []
            residue_code_insertion_list = []
            x_list = []
            y_list = []
            z_list = []
            atom_occupancy_list = []
            atom_bfactor_list = []
            element_list = []

            # lock occupany factor at 1 (instead of: atom.occupancy)
            locked_occupany_factor = 1.00
            max_no_atoms_in_base10 = 99999  # 99,999 for atoms in psf/pdb

            res_no_chain_iter_corrected = []
            res_chain_iteration_corrected_list = []
            residue_id_list = []
            residue_id_adder_fixed_struct_wo_bonds = (
                0  # for example zeolite used as fixed atoms wo bonds
            )
            for i, atom_iter in enumerate(stuct_only_iteration.atoms):
                if i > 0:
                    if (
                        atom_iter.residue.chain == previous_residue_chain
                        and len(atom_iter.bonds) == 0
                    ):
                        residue_id_adder_fixed_struct_wo_bonds += 1

                previous_residue_chain = atom_iter.residue.chain
                residue_id_int = int(
                    unique_residue_data_dict[residue_data_list[i]]
                    + residue_id_adder_fixed_struct_wo_bonds
                )
                res_chain_iteration_corrected_list.append(
                    base10_to_base26_alph(
                        int(residue_id_int / (self.max_residue_no + 1))
                    )[-1:]
                )
                res_id_adder = int(
                    (residue_id_int % self.max_residue_no) % self.max_residue_no
                )
                if int(res_id_adder) == 0:
                    res_no_chain_iter_corrected.append(int(self.max_residue_no))
                else:
                    res_no_chain_iter_corrected.append(res_id_adder)

                residue_id_list.append(residue_id_int)

            # This converts the atom name in the CHARMM psf and pdb files to unique atom names
            [
                unique_individual_atom_names_dict,
                individual_atom_names_list,
                missing_bead_to_atom_name,
            ] = unique_atom_naming(
                stuct_only_iteration,
                residue_id_list,
                residue_names_list,
                bead_to_atom_name_dict=self.bead_to_atom_name_dict,
            )

            if None in [
                unique_individual_atom_names_dict,
                individual_atom_names_list,
                missing_bead_to_atom_name,
            ]:
                self.input_error = True
                print_error_message = (
                    "ERROR: The unique_atom_naming function failed while "
                    "running the charmm_writer function. Ensure the proper inputs are "
                    "in the bead_to_atom_name_dict."
                )

                raise ValueError(print_error_message)

            for coord_iter, atom_coordinates in enumerate(all_atom_coordinates):

                for PDB_residue_count in stuct_only_iteration.residues:
                    segment_id = ""
                    atom_iteration = sorted(
                        PDB_residue_count.atoms, key=lambda atom: atom.number
                    )
                    for atom_iteration_2 in atom_iteration:
                        x, y, z = atom_coordinates[atom_iteration_2.idx]
                        atom_alternate_location_list.append(
                            atom_iteration_2.altloc
                        )
                        residue_code_insertion_list.append(
                            PDB_residue_count.insertion_code[:1]
                        )
                        x_list.append(x)
                        y_list.append(y)
                        z_list.append(z)
                        atom_occupancy_list.append(atom_iteration_2.occupancy)
                        atom_bfactor_list.append(atom_iteration_2.bfactor)
                        element_list.append(
                            Element[atom_iteration_2.atomic_number].upper()
                        )

                for v, atom_iter_1 in enumerate(stuct_only_iteration.atoms):

                    if v + 1 > max_no_atoms_in_base10:
                        atom_number = base10_to_base16_alph_num(v + 1)

                    else:
                        atom_number = v + 1

                    output_write.write(
                        pdb_atom_line_format
                        % (
                            atom_number,
                            individual_atom_names_list[v],
                            atom_alternate_location_list[v],
                            str(residue_names_list[v])[: self.max_resname_char],
                            res_chain_iteration_corrected_list[v],
                            res_no_chain_iter_corrected[v],
                            residue_code_insertion_list[v],
                            x_list[v],
                            y_list[v],
                            z_list[v],
                            locked_occupany_factor,
                            fix_atoms_list[v],
                            segment_id,
                            element_list[v],
                            "",
                        )
                    )

                output_write.write("%-80s\n" % "END")

            output_write.close()

            # **********************************
            # **********************************
            # pdb writer (end)
            # **********************************
            # **********************************
