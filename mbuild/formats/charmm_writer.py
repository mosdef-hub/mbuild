from collections import OrderedDict
from warnings import warn
import os
import numpy as np
from mbuild import Box
from mbuild.utils.conversion import RB_to_CHARMM
from mbuild.utils.sorting import natural_sort
from mbuild.utils.conversion import base10_to_base16_alph_num
from mbuild.utils.conversion import base10_to_base62_alph_num
from mbuild.utils.specific_FF_to_residue import Specific_FF_to_residue
import datetime

from parmed.utils.io import genopen
from parmed.periodic_table import Element
from parmed.charmm.psf import set_molecules
from parmed.utils.six import string_types
from parmed.utils.six.moves import range



def _number_truncated_to_n_digits(num, digits):
    """ Truncates the given number to the specified number of digits """
    if num < 0:
        return int(-(-num % eval('1e%d' % (digits-1))))
    return int(num % eval('1e%d' % digits))

def print_atoms(atom, coords):
    return atom, atom.other_locations, coords[atom.idx]

def _get_bond_types(structure, bonds, sigma_conversion_factor,
        epsilon_conversion_factor):
    unique_bond_types = dict(enumerate(set([(round(bond.type.k*(
                                             sigma_conversion_factor**2/epsilon_conversion_factor),3),
                                             round(bond.type.req/sigma_conversion_factor,3),
                                             tuple(sorted((bond.atom1.type,bond.atom2.type))),
                                             bond.atom1.residue.name, bond.atom2.residue.name
                                             ) for bond in structure.bonds])))

    unique_bond_types = OrderedDict([(y,x+1) for x,y in unique_bond_types.items()])
    bond_types = [unique_bond_types[(round(bond.type.k*(sigma_conversion_factor**2/epsilon_conversion_factor),3),
                                     round(bond.type.req/sigma_conversion_factor,3),
                                     tuple(sorted((bond.atom1.type,bond.atom2.type))),
                                     bond.atom1.residue.name, bond.atom2.residue.name
                                     )] for bond in structure.bonds]

    unique_bond_check_dict = {}
    for i_value_bond, i_key_bond in unique_bond_types.items():
        i_value_duplicated = False
        for j_value_bond, j_key_bond in unique_bond_types.items():
            if j_key_bond > i_key_bond:
                j_value_bond_reorder = (j_value_bond[0], j_value_bond[1],
                                       j_value_bond[2][0], j_value_bond[2][0],
                                        j_value_bond[3],  j_value_bond[4])

                if i_value_bond == j_value_bond_reorder:
                    i_value_duplicated = True
                    if i_value_bond[2][0] > j_value_bond[2][0]:
                        unique_bond_check_dict.update({j_value_bond: len(unique_bond_check_dict) })
                    else:
                        unique_bond_check_dict.update({i_value_bond: len(unique_bond_check_dict) })

            if i_value_duplicated == False:
                unique_bond_check_dict.update({i_value_bond: len(unique_bond_check_dict)})

    unique_bond_types = OrderedDict([(y, x) for y, x in unique_bond_check_dict.items()])

    return bond_types, unique_bond_types

def _get_angle_types(structure, use_urey_bradleys,
        sigma_conversion_factor, epsilon_conversion_factor):
    if use_urey_bradleys:
        charmm_angle_types = []
        for angle in structure.angles:
            ub_k = 0
            ub_req = 0
            for ub in structure.urey_bradleys:
                if (angle.atom1, angle.atom3) == (ub.atom1, ub.atom2):
                    ub_k = ub.type.k
                    ub_req = ub.type.req
            charmm_angle_types.append((round(angle.type.k*(
                sigma_conversion_factor**2/epsilon_conversion_factor),3),
                                       round(angle.type.theteq,3),
                                       round(ub_k/epsilon_conversion_factor, 3),
                                       round(ub_req, 3),
                                       tuple(sorted((angle.atom1.type,angle.atom3.type))),
                                       angle.atom1.residue.name, angle.atom3.residue.name))

        unique_angle_types = dict(enumerate(set(charmm_angle_types)))
        unique_angle_types = OrderedDict([(y,x+1) for x,y in unique_angle_types.items()])
        angle_types = [unique_angle_types[ub_info] for ub_info in charmm_angle_types]

    else:
        unique_angle_types = dict(enumerate(set([(round(angle.type.k*(
            sigma_conversion_factor**2/epsilon_conversion_factor),3),
                                                  round(angle.type.theteq,3),
                                                  angle.atom2.type,
                                                  tuple(sorted((angle.atom1.type,angle.atom3.type))),
                                                  angle.atom1.residue.name, angle.atom2.residue.name,
                                                  angle.atom3.residue.name
                                                  ) for angle in structure.angles])))
        unique_angle_types = OrderedDict([(y,x+1) for x,y in unique_angle_types.items()])
        angle_types = [unique_angle_types[(round(angle.type.k*(
            sigma_conversion_factor**2/epsilon_conversion_factor),3),
                                           round(angle.type.theteq,3),
                                           angle.atom2.type,
                                           tuple(sorted((angle.atom1.type,angle.atom3.type))) ,
                                           angle.atom1.residue.name, angle.atom2.residue.name,
                                           angle.atom3.residue.name
                                           )] for angle in structure.angles]

    unique_angle_check_dict = {}
    for i_value_ang, i_key_ang in unique_angle_types.items():
        i_value_duplicated = False
        for j_value_ang, j_key_ang in unique_angle_types.items():
            if j_key_ang > i_key_ang:
                j_value_ang_reorder = (j_value_ang[0], j_value_ang[1],
                                       j_value_ang[2], j_value_ang[3][0], j_value_ang[3][1],
                                         j_value_ang[4], j_value_ang[5], j_value_ang[6])

                if i_value_ang == j_value_ang_reorder:
                    i_value_duplicated = True
                    if i_value_ang[2] > j_value_ang[2]:
                        unique_angle_check_dict.update({j_value_ang: len(unique_angle_check_dict) })
                    else:
                        unique_angle_check_dict.update({i_value_ang: len(unique_angle_check_dict) })

            if i_value_duplicated == False:
                unique_angle_check_dict.update({i_value_ang: len(unique_angle_check_dict)})

    unique_angle_types = OrderedDict([(y, x) for y, x in unique_angle_check_dict.items()])

    return angle_types, unique_angle_types

def _get_dihedral_types(structure, use_rb_torsions, use_dihedrals,
         epsilon_conversion_factor):
    lj_unit = 1 / epsilon_conversion_factor
    if use_rb_torsions:
        unique_dihedral_types = dict(enumerate(set([(round(dihedral.type.c0*lj_unit,3),
                                                     round(dihedral.type.c1*lj_unit,3),
                                                     round(dihedral.type.c2*lj_unit,3),
                                                     round(dihedral.type.c3*lj_unit,3),
                                                     round(dihedral.type.c4*lj_unit,3),
                                                     round(dihedral.type.c5*lj_unit,3),
                                                     round(dihedral.type.scee,1),
                                                     round(dihedral.type.scnb,1),
                                                     dihedral.atom1.type, dihedral.atom2.type,
                                                     dihedral.atom3.type, dihedral.atom4.type,
                                                     dihedral.atom1.residue.name, dihedral.atom2.residue.name,
                                                     dihedral.atom3.residue.name, dihedral.atom4.residue.name
                                                     ) for dihedral in structure.rb_torsions])))



        unique_dihedral_types = OrderedDict([(y,x+1) for x,y in unique_dihedral_types.items()])

        dihedral_types = [unique_dihedral_types[(round(dihedral.type.c0*lj_unit,3),
                                                 round(dihedral.type.c1*lj_unit,3),
                                                 round(dihedral.type.c2*lj_unit,3),
                                                 round(dihedral.type.c3*lj_unit,3),
                                                 round(dihedral.type.c4*lj_unit,3),
                                                 round(dihedral.type.c5*lj_unit,3),
                                                 round(dihedral.type.scee,1),
                                                 round(dihedral.type.scnb,1),
                                                 dihedral.atom1.type, dihedral.atom2.type,
                                                 dihedral.atom3.type, dihedral.atom4.type,
                                                 dihedral.atom1.residue.name, dihedral.atom2.residue.name,
                                                 dihedral.atom3.residue.name, dihedral.atom4.residue.name
                                                 )] for dihedral in structure.rb_torsions]

    elif use_dihedrals:
        charmm_dihedrals = []
        structure.join_dihedrals()
        for dihedral in structure.dihedrals:
            if not dihedral.improper:
                weight = 1 / len(dihedral.type)
                for dih_type in dihedral.type:
                    charmm_dihedrals.append((round(dih_type.phi_k*lj_unit,3),
                                             int(round(dih_type.per,0)),
                                             int(round(dih_type.phase,0)),
                                             round(weight, 4),
                                             round(dih_type.scee,1),
                                             round(dih_type.scnb,1),
                                             dihedral.atom1.type, dihedral.atom2.type,
                                             dihedral.atom3.type, dihedral.atom4.type,
                                             dihedral.atom1.residue.name, dihedral.atom2.residue.name,
                                             dihedral.atom3.residue.name, dihedral.atom4.residue.name
                                             ))

        unique_dihedral_types = dict(enumerate(set(charmm_dihedrals)))
        unique_dihedral_types = OrderedDict([(y,x+1) for x,y in unique_dihedral_types.items()])
        dihedral_types = [unique_dihedral_types[dihedral_info] for dihedral_info in charmm_dihedrals]

    unique_dihedral_check_dict = OrderedDict()
    for i_value_dihed, i_key_dihed in unique_dihedral_types.items():
        i_value_duplicated = False
        for j_value_dihed, j_key_dihed in unique_dihedral_types.items():
            if j_key_dihed > i_key_dihed:
                j_value_dihed_reorder = (j_value_dihed[0], j_value_dihed[1], j_value_dihed[2], j_value_dihed[3],
                                         j_value_dihed[4], j_value_dihed[5], j_value_dihed[6], j_value_dihed[7],
                                         j_value_dihed[11], j_value_dihed[10], j_value_dihed[9], j_value_dihed[8],
                                         j_value_dihed[15], j_value_dihed[14], j_value_dihed[13], j_value_dihed[12])

                if i_value_dihed == j_value_dihed_reorder:
                    i_value_duplicated = True
                    if i_value_dihed[8] > j_value_dihed[8]:
                        unique_dihedral_check_dict.update({j_value_dihed: len(unique_dihedral_check_dict)+1})
                    else:
                        unique_dihedral_check_dict.update({i_value_dihed: len(unique_dihedral_check_dict)+1})

        if i_value_duplicated == False:
            unique_dihedral_check_dict.update({i_value_dihed: len(unique_dihedral_check_dict)+1 })

    unique_dihedral_types = OrderedDict([(y, x) for y, x in unique_dihedral_check_dict.items()])

    return dihedral_types, unique_dihedral_types

def _get_impropers(structure, epsilon_conversion_factor):
    lj_unit = 1 / epsilon_conversion_factor
    unique_improper_types = dict(enumerate(set([(round(improper.type.psi_k*lj_unit,3),
                                                 round(improper.type.psi_eq,3),
                                                 improper.atom1.type, improper.atom2.type,
                                                 improper.atom3.type, improper.atom4.type,
                                                 improper.atom1.residue.name, improper.atom2.residue.name,
                                                 improper.atom3.residue.name, improper.atom4.residue.name)
                                                for improper in structure.impropers])))
    unique_improper_types = OrderedDict([(y,x+1) for x,y in unique_improper_types.items()])
    improper_types = [unique_improper_types[(round(improper.type.psi_k*lj_unit,3),
                                             round(improper.type.psi_eq,3),
                                             improper.atom1.type, improper.atom2.type,
                                             improper.atom3.type, improper.atom4.type,
                                             improper.atom1.residue.name, improper.atom2.residue.name,
                                             improper.atom3.residue.name, improper.atom4.residue.name)]
                      for improper in structure.impropers]

    # If impropers are added to GOMC, add the atom sorter for the unique combinations here

    return improper_types, unique_improper_types


def unique_atom_naming(structure, residue_ID_list, residue_names_list, Bead_to_atom_name_dict=None):

    """ Outputs
        unique_Individual_atom_names_dict : dictionary
            All the unique atom names compiled into a dictionary.
        Individual_atom_names_List : list, in sequential  order
            The atom names for every atom in the system
        Missing_Bead_to_atom_name  : list, in sequential  order
            The bead names of any atoms beads that did not have a name specificed to them
            via the Bead_to_atom_name_dict

    Parameters
    ----------
    structure : compound object
    residue_ID_list : list, in sequential  order
            The residue ID for every atom in the system
    residue_names_list  : list, in sequential  order
        The atom names for every atom in the system
    Bead_to_atom_name_dict: dictionary ; optional, default =None
        For all atom names/elements/beads with 2 or less digits, this converts
        the atom name in the GOMC psf and pdb files to a unique atom name,
        provided they do not exceed 3844 atoms (62^2) of the same name/element/bead
        per residue. For all atom names/elements/beads with 3 digits, this converts
        the atom name in the GOMC psf and pdb files to a unique atom name,
        provided they do not exceed 62 of the same name/element pre residue.
        Example dictionary: {'_CH3':'C', '_CH2':'C', '_CH':'C', '_HC':'C'}
    """
    unique_Individual_atom_names_dict = {}
    Individual_atom_names_List = []
    Missing_Bead_to_atom_name = []
    for i, atom in enumerate(structure.atoms):
        interate_thru_names = True
        j = 0
        while interate_thru_names == True:
            j = j + 1
            if str(atom.name)[:1] == '_':
                if Bead_to_atom_name_dict != None and (str(atom.name) in Bead_to_atom_name_dict) == True:
                    if len(Bead_to_atom_name_dict[str(atom.name)]) > 2:
                        text_to_write = ('ERROR: only enter atom names that have 2 or less digits' +
                                         ' in the Bead to atom naming dictionary (Bead_to_atom_name_dict) ')
                        return print(warn(text_to_write))
                    else:
                        atom_name_value = Bead_to_atom_name_dict[str(atom.name)]
                        No_digits_atom_name = 2
                else:
                    Missing_Bead_to_atom_name.append(1)
                    atom_name_value = 'BD'
                    No_digits_atom_name = 2
            elif len(str(atom.name)) > 2:
                if len(str(atom.name)) == 3:
                    No_digits_atom_name = 1
                else:
                    text_to_write = ('ERROR: atom numbering will not work propery at' +
                                     ' the element has more than 4 charaters')
                    return print(warn(text_to_write))
            else:
                No_digits_atom_name = 2
                atom_name_value = atom.name
            atom_name_iteration = str(atom_name_value) + str(base10_to_base62_alph_num(j))
            atom_ResNo_ResName_AtomName_iteration = str(residue_ID_list[i]) + '_' \
                                                    + str(residue_names_list[i]) + '_' + atom_name_iteration

            if unique_Individual_atom_names_dict.get(str(atom_ResNo_ResName_AtomName_iteration)) is None:
                unique_Individual_atom_names_dict.update({atom_ResNo_ResName_AtomName_iteration: i + 1})
                interate_thru_names = False
                Individual_atom_names_List.append(
                    str(atom_name_value) + str(str(base10_to_base62_alph_num(j))[-No_digits_atom_name:]))

    if sum(Missing_Bead_to_atom_name) > 0:
        warn("NOTE: All bead names were not found in the Bead to atom naming dictionary (Bead_to_atom_name_dict) ")

    return unique_Individual_atom_names_dict, Individual_atom_names_List, Missing_Bead_to_atom_name


# Currently the NBFIX is disabled as since only the OPLS and TRAPPE force fields are currently supported

def charmm_psf_psb_FF(structure_0, filename_0, structure_1 = None, filename_1= None,
                      non_bonded_type='LJ', forcefield_files=None, forcefield_names = None,  residues=None,
                      detect_forcefield_style=True, fix_res_bonds_angles = None, Bead_to_atom_name_dict=None,
                      fix_residue=None, fix_residue_in_box=None,  coordinates=None, GOMC_FF_filename= None,
                      reorder_res_in_pdb_psf =False, box_0 = None, box_1 = None  ,**kwargs):

    """Output a GOMC data file.

    Outputs a GOMC data file The units are as follows
        * Mw = g/mol
        * Harmonic bonds : Kb = kcal/mol, b0 = Angstroms
        * Harmonic angles : Ktheta = kcal/mole/rad**2 , Theta0 = degrees
        * Dihedral angles: Ktheta = kcal/mole, n = interger (unitless), delta = degrees
        * nonbonded : epsilon = kcal/mol, Rmin = Angstroms, n = interger (unitless)
        Note: units are the same at the LAMMPS real units.  The atom style
        is the same as the lammps 'full' atom style format.

    Parameters
    ----------
    structure_0 : compound object
    filename_0 : str
        Path of the output file for structure_0
    structure_1 : compound object number 2, optional
        (Ex: for GCMC or GEMC simulations which have mulitiple simulation boxes)
    filename_1 : str , optional
        Path of the output file for structure_1 (Ex: for GCMC or GEMC simulations
        which have mulitiple simulation boxes)
    non_bonded_type : str. optional, default = 'LJ' (i.e., Lennard-Jones )
            Specify the type of non-bonded potential for the GOMC force field files.
            Note: Currently, on the 'LJ' potential is supported.
    residues : str of list of str
        Labels of residues in the Compound. Residues are assigned by
        checking against Compound.name.
        Note: to write the GOMC force field files and the psf files the
        residues and forcefield_files/forcefiled_names must be provided in a list, which
        is in sequential order of each other. (Example:
        residues = [WAT, OIL] and forcefield_files/forcefiled_names = [WAT_FF_file, OIL_FF_file]
    forcefield_files : str or dictionary, default = None
        Apply a forcefield to the output file using a forcefield provided
        by the `foyer` package.
        Note: to write the GOMC force field files and the psf files the
        residues and forcefield_files must be provided in a str or
        dictionary.  If a dictionary is provided all residues must
        be specified to a force field.
            Ex dict: {'Water' : 'oplsaa.xml', 'OCT': 'trappe-ua.xml'}
            Ex str: 'trappe-ua.xml'
        Note: the file path name must also be specified.
    forcefield_names : str or dictionary, default = None
        Apply a forcefield to the output file using a forcefield provided
        by the `foyer` package.  These forcefields are called by name and
        stored in the prebuilt foyer package.
        Note: to write the GOMC force field files and the psf files the
        residues and forcefield_names must be provided in a str or
        dictionary.  If a dictionary is provided all residues must
        be specified to a force field.
            Ex dict: {'Water' : 'oplsaa', 'OCT': 'trappe-ua'}
            Ex str: 'trappe-ua'
        Note: the file path name must also be specified.
    detect_forcefield_style: boolean
        If True, format lammpsdata parameters based on the contents of
        the parmed structure_0
    fix_res_bonds_angles: list, default = None
        When list of residues is provided the the residue will have
        bonds and angles fixed.  This is only for the GOMC force field
        writer.
        WARNING: Currently if the residue is similar to another residue,
        this will not work as distinguishing between the 2 residues is
        currently under construction.
        Example 1: If you have a n-pentane and n-octane, you will not
        be able to fix the bonds and angles on one and not the other
        Example 2: if you have a water and n-pentane, you are able to
        fix the waters and/or n-pentanes bonds and angles.
    Bead_to_atom_name_dict: dict, optional, default =None
        For all atom names/elements/beads with 2 or less digits, this converts
        the atom name in the GOMC psf and pdb files to a unique atom name,
        provided they do not exceed 3844 atoms (62^2) of the same name/element/bead
        per residue. For all atom names/elements/beads with 3 digits, this converts
        the atom name in the GOMC psf and pdb files to a unique atom name,
        provided they do not exceed 62 of the same name/element pre residue.
        Example dictionary: {'_CH3':'C', '_CH2':'C', '_CH':'C', '_HC':'C'}
    fix_residue: list  or None, default = None
        Changes occcur in the GOMC pdb file only.
        When residues are listed here, all the atoms in the residue are
        fixed and can not move via setting the Beta values in the PDB
        file to 1.00.
        If neither fix_residue or fix_residue_in_box lists a
        residue or both equal None, then the Beta values for all the atoms
        in the residue are free to move in the simulation and Beta values
        in the PDB file is set to 0.00
    fix_residue_in_box: list  or None, default = None
        Changes occcur in the GOMC pdb file only.
        When residues are listed here, all the atoms in the residue become
        can move within the box but cannot be transferred between boxes
        via setting the Beta values in the PDB file to 2.00.
        If neither fix_residue or fix_residue_in_box lists a
        residue or both equal None, then the Beta values for all the atoms
        in the residue are free to move in the simulation and Beta values
        in the PDB file is set to 0.00
    GOMC_FF_filename ; str, default =None
        If a sting, it will write the GOMC force field files for the
        structures
    reorder_res_in_pdb_psf ; bool, default =False
        If False, the order of of the atoms in the pdb file is kept in
        its original order, as in the Compound sent to the writer.
        If True, the order of the atoms is reordered based on their
        residue names in the 'residues' list that was entered.
    box_0 ; list of 3 positive float values or the dimensions [x, y ,z]
        for structure_0 in nanometers (nm)
        This is to add/override or change the structures dimenstions. Ex: [1,2,3]
    box_1 ; list of 3 positive float values or the dimensions [x, y ,z]
        for structure_1 in nanometers (nm)
        This is to add/override or change the structures dimenstions. Ex: [1,2,3]
    Notes
    -----
    Impropers and NBFIX are not currenly supported
    Currently the NBFIX is disabled as since only the OPLS and TRAPPE force fields are supported
    OPLS and CHARMM forcefield styles are supported, AMBER forcefield styles are NOT
    Impropers and Urey-Bradleys are not supported for GOMC

    The atom typing is currently numbering.  Therefore, if you are utilizing
    multiple/boxes (i.e., GCMC or GEMC), you need to ensure that all the same molecules
    are in each system or there are zero molecules in a systems.  This will ensure that
    they atom typing all matches and works properly with the FF, pdb, and psf files.
    The atom typing will be changed to unique naming to avoid this issue in the future.

    Unique atom names are provided if the system do not exceed 3844 atoms (62^2) of the same
    name/bead per residue. For all atom names/elements with 3 or less digits, this converts
    the atom name in the GOMC psf and pdb files to a unique atom name,  provided they do
    not exceed 62 of the same name/element pre residue.

    Generating an empty pdb/psf:
        Single System: Enter residues = [], but the accompanying structure (structure_0)
        must be an empty mb.Compound structure. However, when doing this, the forcefield_names or
        forcefield_files must be supplied, or it will provide an error
        (i.e., forcefield_files and forcefield_names are both not equal to None)
        Dual System: Enter an empty mb.Compound structure for either structure_0 or
        structure_1.

    """


    if not isinstance(residues, list):
        return warn('Error: Please enter the residues (residues) in a list format')

    date_time = datetime.datetime.today()

    if residues is None:
        return warn('Error: Please enter the residues (residues)  list')
    if not isinstance(filename_0, str):
        return warn('Error: Please enter the filename_0 as a string')
    if filename_1 != None and not isinstance(filename_1, str):
        return warn('Error: Please enter the filename_1 as a string')

    if structure_1 is None and box_1 != None:
        warn('Error: box_1 is set to a value but there is not a structure 1 to use it on.')


    if GOMC_FF_filename != None:
        if not isinstance(GOMC_FF_filename, str) :
            return warn('Error: Please enter GOMC force field name (GOMC_FF_filename) as a string')
        if isinstance(GOMC_FF_filename, str) :
            extension_FF_name = os.path.splitext(GOMC_FF_filename)[-1]
            if extension_FF_name == '.inp':
                print('')
            else:
                GOMC_FF_filename = GOMC_FF_filename + '.inp'

    if forcefield_names is None and forcefield_files is None:
        return warn('Please enter either the forcefield_files or forcefield_names, neither were provided')

    if forcefield_names != None and forcefield_files != None:
        return warn('Please enter either the forcefield_files or forcefield_names, not both')

    if forcefield_files != None and  forcefield_names is None:

        print('write_gomcdata: forcefield_files = '+str(forcefield_files) +', ' +'residues = '+str(residues) )

        if forcefield_files != None and not isinstance(forcefield_files, dict) and not isinstance(forcefield_files, str):
            return warn('The force field file (forcefield_files) is not a string or a dictionary with'+
                        ' all the residues specified to a force field.'+
                          "-> String Ex: 'path/trappe-ua.xml'  ."+
                         "Otherise provided a dictionary with all the residues specified to a force field "+
                         "->Dictionary Ex: {'Water' : 'path/oplsaa.xml', 'OCT': 'path/trappe-ua.xml'}, "+
                         "Note: the file path must be specified the force field file")





        if isinstance(forcefield_files, list) == True:
            return warn('Error: Please enter the forcefield_files (forcefield_files) as a single string' +
                        '  or in a dictionary format (i.e., residue to for each FF)')

        if isinstance(forcefield_files, str)== True :
            FF_name = forcefield_files
            forcefield_files = {}
            for i in range(0, len(residues)):
                forcefield_files.update({residues[i] : FF_name})
            print('FF forcefield_files = '+str(forcefield_files))

    elif forcefield_names != None and  forcefield_files is None:
        print('write_gomcdata: forcefield_names = ' + str(forcefield_names) + ', ' + 'residues = ' + str(residues))

        if forcefield_names != None and not isinstance(forcefield_names, dict) and not isinstance(forcefield_names,str):
            return warn('The force field names (forcefield_names) is not a string or a dictionary with' +
                        ' all the residues specified to a force field.' +
                        "-> String Ex: 'trappe-ua' or 'oplsaa'  ." +
                        "Otherise provided a dictionary with all the residues specified to a force field " +
                        "->Dictionary Ex: {'Water' : 'oplsaa', 'OCT': 'trappe-ua'}, " +
                        "Note: the file path must be specified the force field file")

        if isinstance(forcefield_names, list) == True:
            return warn('Error: Please enter the forcefield_names (forcefield_names) as a single string' +
                        '  or in a dictionary format (i.e., residue to for each FF)')

        if isinstance(forcefield_names, str) == True:
            FF_name = forcefield_names
            forcefield_names = {}
            for i in range(0, len(residues)):
                forcefield_names.update({residues[i]: FF_name})
            print('FF forcefield_names = ' + str(forcefield_names))

    if residues != None and not isinstance(residues, list):
        return warn('Error: Please enter the residues (residues) in a list format')




    if fix_res_bonds_angles != None and not isinstance(fix_res_bonds_angles, list):
        return warn('Error: Please enter the residues that have fixed angles' +
                    ' and bonds (fix_res_bonds_angles) in a list format')

    if isinstance(fix_res_bonds_angles, list):
        for q in range(0,len(fix_res_bonds_angles)):
            if fix_res_bonds_angles[q] not in residues:
                return warn('Error: Please ensure that all the residue names in the fix_res_bonds_angles '+
                            'list are also in the residues list.')
            else:
                print('INFORMATION: The following residues will have fixed bonds'
                      + ' and angles: fix_res_bonds_angles = ' +str(fix_res_bonds_angles))


    if fix_residue != None and not isinstance(fix_residue, list):
        return warn('Error: Please enter the fix_residue in a list format')

    if isinstance(fix_residue, list):
        for q in range(0,len(fix_residue)):
            if fix_residue[q] not in residues:
                return warn('Error: Please ensure that all the residue names in the fix_residue '+
                            'list are also in the residues list.')

    if fix_residue_in_box != None and not isinstance(fix_residue_in_box, list):
        return warn('Error: Please enter the fix_residue_in_box in a list format')

    if isinstance(fix_residue_in_box, list):
        for q in range(0,len(fix_residue_in_box)):
            if fix_residue_in_box[q] not in residues:
                return warn('Error: Please ensure that all the residue names in the fix_residue_in_box '+
                            'list are also in the residues list.')

    if Bead_to_atom_name_dict != None and not isinstance(Bead_to_atom_name_dict, dict):
        return warn('Error: Please enter the a bead type to atom in the' +
                    '  dictionary (Bead_to_atom_name_dict) '
                    + 'so GOMC can properly evaluate the unique atom names')


    if box_0 !=None :
        box_length = len(box_0)
        if box_length != 3:
             return warn('Please enter all 3 values for the box_0 dimensions.')
        for box_iter in range(0, len(box_0)):
            if isinstance(box_0[ box_iter], str)==True:
                return warn('Please enter all positive or 0 values for the box_0 dimensions.')
            if box_0[ box_iter] < 0:
                return warn('Please enter all positive or 0 values for the box_0 dimensions.')

    if box_1 !=None :
        box_length = len(box_1)
        if box_length != 3:
             return warn('Please enter all 3 values for the box_1 dimensions.')
        for box_iter in range(0, len(box_1)):
            if isinstance(box_1[ box_iter], str)==True:
                return warn('Please enter all positive or 0 values for the box_1 dimensions.')
            if box_1[ box_iter] < 0:
                return warn('Please enter all positive or 0 values for the box_1 dimensions.')


    print("******************************")
    print("")

    if structure_1 != None:
        boxes_for_simulation = 2
    else:
        boxes_for_simulation = 1

    #write the Force fields
    combined_1_4_LJ_dict_per_residue = {}
    combined_1_4_Coul_dict_per_residue = {}
    if structure_1 != None:

        print('GOMC FF writing each residues FF as a group for structure_0')
        structure_0_FF, \
        coulomb14scaler_dict_0, \
        LJ14scaler_dict_0,\
        residues_applied_list_0 = Specific_FF_to_residue(structure_0,
                                                         forcefield_files=forcefield_files,
                                                         forcefield_names=forcefield_names,
                                                         residues=residues,
                                                         reorder_res_in_pdb_psf=reorder_res_in_pdb_psf,
                                                         box = box_0,
                                                         boxes_for_simulation = boxes_for_simulation)
        print('GOMC FF writing each residues FF as a group for  structure_1')
        structure_1_FF, \
        coulomb14scaler_dict_1, \
        LJ14scaler_dict_1, \
        residues_applied_list_1 = Specific_FF_to_residue(structure_1,
                                                         forcefield_files=forcefield_files,
                                                         forcefield_names=forcefield_names,
                                                         residues=residues,
                                                         reorder_res_in_pdb_psf=reorder_res_in_pdb_psf,
                                                         box = box_1,
                                                         boxes_for_simulation = boxes_for_simulation)
        structure_0_and_1_FF =structure_0_FF + structure_1_FF
        combined_1_4_LJ_dict_per_residue.update(coulomb14scaler_dict_0)
        combined_1_4_LJ_dict_per_residue.update(coulomb14scaler_dict_1)
        combined_1_4_Coul_dict_per_residue.update(coulomb14scaler_dict_0)
        combined_1_4_Coul_dict_per_residue.update(coulomb14scaler_dict_1)

        residues_applied_list_0_and_1 = residues_applied_list_0
        for res_iter in range(0,len(residues_applied_list_1 )):
            if residues_applied_list_1[res_iter] not in residues_applied_list_0:
                residues_applied_list_0_and_1.append(residues_applied_list_1[res_iter])

        for res_iter_1 in range(0, len(residues_applied_list_0_and_1)):
            if residues_applied_list_0_and_1[res_iter_1] not in  residues:
                return warn("All the residues were not used from the forcefield_names or forcefield_files "+ \
                            "string or dictionary.  There may be residues below other specified residues "+ \
                            "in the mbuild.Compound hierarchy.  If so, the residues acquire the residue's "+ \
                            "force fields, which is at the top of the hierarchy.  Alternatively, "+ \
                            "residues that are not in the structure may have been specified.")

        total_charge = sum([atom.charge for atom in structure_0_FF])
        if round(total_charge, 4) != 0.0:
            warn('System is not charge neutral for structure_0. Total charge is {}.'.format(total_charge))

        total_charge = sum([atom.charge for atom in structure_1_FF])
        if round(total_charge, 4) != 0.0:
            warn('System is not charge neutral for structure_1. Total charge is {}.'.format(total_charge))

        total_charge = sum([atom.charge for atom in structure_0_and_1_FF])
        if round(total_charge, 4) != 0.0:
            warn('System is not charge neutral for structure__0_and_1. Total charge is {}.'.format(total_charge))

    else:


        print('GOMC FF writing each residues FF as a group for structure_0')
        structure_0_FF, \
        coulomb14scaler_dict_0, \
        LJ14scaler_dict_0, \
        residues_applied_list_0 = Specific_FF_to_residue(structure_0,
                                                         forcefield_files=forcefield_files,
                                                         forcefield_names=forcefield_names,
                                                         residues=residues,
                                                         reorder_res_in_pdb_psf=reorder_res_in_pdb_psf,
                                                         box=box_0,
                                                         boxes_for_simulation = boxes_for_simulation)
        combined_1_4_LJ_dict_per_residue.update(coulomb14scaler_dict_0)
        combined_1_4_Coul_dict_per_residue.update(coulomb14scaler_dict_0)

        for res_iter_1 in range(0, len(residues_applied_list_0)):
            if residues_applied_list_0[res_iter_1] not in residues:
                return warn("All the residues were not used from the forcefield_names or forcefield_files " + \
                            "string or dictionary.  There may be residues below other specified residues " + \
                            "in the mbuild.Compound hierarchy.  If so, the residues acquire the residue's " + \
                            "force fields, which is at the top of the hierarchy.  Alternatively, " + \
                            "residues that are not in the structure may have been specified.")

        total_charge = sum([atom.charge for atom in structure_0_FF])
        if round(total_charge, 4) != 0.0:
            warn('System is not charge neutral for structure_0. Total charge is {}.'.format(total_charge))


    unique_residue_data_dict = {}
    unique_residue_data_list = []
    residue_data_name_list = []
    if structure_1 != None:
        residue_iterate = 0
        for m, residue in enumerate(structure_0_FF.residues):
            residue_iterate = residue_iterate + 1
            unique_residue_data_list.append(str(structure_0_FF.residues[m]))
            unique_residue_data_dict.update({unique_residue_data_list[m]: m + 1})
            residue_data_name_list.append(structure_0_FF.residues[m].name)

        for m, residue in enumerate(structure_1_FF.residues):
            unique_residue_data_list.append(str(structure_1_FF.residues[m]))
            unique_residue_data_dict.update({unique_residue_data_list[m]: m + 1 + residue_iterate})
            residue_data_name_list.append(structure_1_FF.residues[m].name)


    else:
        for m, residue in enumerate(structure_0_FF.residues):
            unique_residue_data_list.append(str(structure_0_FF.residues[m]))
            unique_residue_data_dict.update({unique_residue_data_list[m]: m + 1})
            residue_data_name_list.append(structure_0_FF.residues[m].name)




    for n in range(0, len(residue_data_name_list)):
        if residue_data_name_list[n] not in residues:
            print('residue_data_name_list = '+str(residue_data_name_list))
            return warn('Error: Please specifiy all residues (residues) in a list')

    if forcefield_files != None:
        print('forcefield type from compound = '+str( forcefield_files))
    elif forcefield_names != None:
        print('forcefield type from compound = ' + str(forcefield_names))
    print('coulomb14scale from compound = ' + str(combined_1_4_Coul_dict_per_residue))
    print('lj14scale from compound = ' + str(combined_1_4_LJ_dict_per_residue))



    """
    Note:
    -----
    unique_types : a sorted list of unique atomtypes for all atoms in the structure_0_FF.
        Defined by:
            atomtype : atom.type
    unique_bond_types: an enumarated OrderedDict of unique bond types for all bonds in the structure_0_FF.
        Defined by bond parameters and component atomtypes, in order:
            k : bond.type.k
            req : bond.type.req
            atomtypes : sorted((bond.atom1.type, bond.atom2.type))
    unique_angle_types: an enumerated OrderedDict of unique angle types for all angles in the structure_0_FF.
        Defined by angle parameters and component atomtypes, in order:
            k : angle.type.k
            theteq : angle.type.theteq
            vertex atomtype: angle.atom2.type
            atomtypes: sorted((bond.atom1.type, bond.atom3.type))
    unique_dihedral_types: an enumerated OrderedDict of unique dihedrals type for all dihedrals in the structure_0_FF.
        Defined by dihedral parameters and component atomtypes, in order:
            c0 : dihedral.type.c0
            c1 : dihedral.type.c1
            c2 : dihedral.type.c2
            c3 : dihedral.type.c3
            c4 : dihedral.type.c4
            c5 : dihedral.type.c5
            scee : dihedral.type.scee
            scnb : dihedral.type.scnb
            atomtype 1 : dihedral.atom1.type
            atomtype 2 : dihedral.atom2.type
            atomtype 3 : dihedral.atom3.type
            atomtype 4 : dihedral.atom4.type
    """
    # lock the atom_style and unit_style for GOMC. Can be inserted into variables once more functionality is built in
    atom_style = 'full'
    unit_style = 'real'
    # functional form type default.  Can be inserted into variables once more functionality is built in
    use_rb_torsions = True
    use_dihedrals = False
    use_urey_bradleys = False

    # Convert coordinates to LJ units
    if unit_style == 'real':
        sigma_conversion_factor = 1
        epsilon_conversion_factor = 1
        mass_conversion_factor = 1
    else:
        return print("unit_style is not real and thus not available to this writer")


    if structure_1 !=None:
        types = [atom.type+'_'+str(atom.residue.name) for atom in structure_0_and_1_FF.atoms]

    else:
        types = [atom.type + '_' + str(atom.residue.name) for atom in structure_0_FF.atoms]

    unique_types = list(set(types))
    unique_types.sort(key=natural_sort)

    print('unique_types = '+str(unique_types))

    #***************
    # edited fix for residues (start)
    # ***************
    if structure_1 !=None:
        masses = np.array([atom.mass for atom in structure_0_and_1_FF.atoms]) / mass_conversion_factor
        mass_dict = dict([(unique_types.index(atom_type) + 1, mass) for atom_type, mass in zip(types, masses)])

    else:
        masses = np.array([atom.mass for atom in structure_0_FF.atoms]) / mass_conversion_factor
        mass_dict = dict([(unique_types.index(atom_type) + 1, mass) for atom_type, mass in zip(types, masses)])


    # added an index so the atom types can be converted to numbers as the type name is to long for insertion into
    # the pdb and psf files
    atom_types_to_index_value_dict = dict(
        [(unique_types[unique_types.index(atom_type)], unique_types.index(atom_type) + 1) for atom_type, masses in
         zip(types, masses)])

    box_0 = Box(lengths=np.array([0.1 * val for val in structure_0_FF.box[0:3]]),angles=structure_0_FF.box[3:6])
    #Divide by conversion factor
    box_0.maxs /= sigma_conversion_factor

    #Internally use nm
    if structure_1 != None:
        box_1 = Box(lengths=np.array([0.1 * val for val in structure_1_FF.box[0:3]]),
                  angles=structure_1_FF.box[3:6])
        # Divide by conversion factor
        box_1.maxs /= sigma_conversion_factor

    if structure_1 != None:
        structure_selection = structure_0_and_1_FF
    else:
        structure_selection = structure_0_FF


    # Syntax which can change based on the functional form
    # Infer functional form based on the properties of the structure_0_and_1_FF or structure_0_FF
    if detect_forcefield_style:
        # Check angles
        if len(structure_selection.urey_bradleys) > 0 :
            print("Urey bradley terms detected")
            data.write("Urey bradley terms detected, will use angle_style charmm")
            data.write("Error GOMC does no support the Urey bradley terms")
            use_urey_bradleys = True
        else:
            print("No urey bradley terms detected, will use angle_style harmonic")
            use_urey_bradleys = False

        # Check dihedrals
        if len(structure_selection.rb_torsions) > 0:
            print("will use CHARMM_torsions  =  K0 + K1 * (1 + Cos[n1*(t) - (d1)] ) + "+
                  "K2 * (1 + Cos[n2*(t) - (d2)] ) + K3 * (1 + Cos[n3*(t) - (d3)] ) + "  +
                  "K4 * (1 + Cos[n4*(t) - (d4)] ) + K5 * (1 + Cos[n5*(t) - (d5)] ) ")
            use_rb_torsions = True


        else:
            use_rb_torsions = False
        if len(structure_selection.dihedrals) > 0:
            print("Charmm dihedrals detected, will use dihedral_style charmm")
            # this will need tested with a standard charmm input format before releasing it
            use_dihedrals = True
            return print("Charmm dihedrals not yet supported ")
        else:
            use_dihedrals = False
    if use_rb_torsions and use_dihedrals:
        raise ValueError("Multiple dihedral styles detected, check your "
                         "Forcefield XML and structure_selection")

    # Check impropers
    for dihedral in structure_selection.dihedrals:
        if dihedral.improper:
            raise ValueError("Amber-style impropers are currently not supported")

    bonds = [[bond.atom1.idx+1, bond.atom2.idx+1] for bond in structure_selection.bonds]
    angles = [[angle.atom1.idx+1,
               angle.atom2.idx+1,
               angle.atom3.idx+1] for angle in structure_selection.angles]
    if use_rb_torsions:
        dihedrals = [[dihedral.atom1.idx+1,
                      dihedral.atom2.idx+1,
                      dihedral.atom3.idx+1,
                      dihedral.atom4.idx+1] for dihedral in structure_selection.rb_torsions]
    elif use_dihedrals:
        dihedrals = [[dihedral.atom1.idx+1,
                      dihedral.atom2.idx+1,
                      dihedral.atom3.idx+1,
                      dihedral.atom4.idx+1] for dihedral in structure_selection.dihedrals]
    else:
        dihedrals = []
    impropers = [[improper.atom1.idx+1,
                  improper.atom2.idx+1,
                  improper.atom3.idx+1,
                  improper.atom4.idx+1] for improper in structure_selection.impropers]


    if bonds :
        if len(structure_selection.bond_types) == 0:
            bond_types = np.ones(len(bonds),dtype=int)
        else:
            bond_types, unique_bond_types = _get_bond_types(structure_selection,
                    bonds, sigma_conversion_factor,
                    epsilon_conversion_factor)


    if angles:
        angle_types, unique_angle_types = _get_angle_types(structure_selection,
                use_urey_bradleys, sigma_conversion_factor,
                epsilon_conversion_factor)

    if dihedrals:
        dihedral_types, unique_dihedral_types = _get_dihedral_types(
                structure_selection, use_rb_torsions, use_dihedrals,
                epsilon_conversion_factor)



    if impropers:
        improper_types, unique_improper_types = _get_impropers(structure_selection,
                epsilon_conversion_factor)

    if GOMC_FF_filename is None:
        print('GOMC FF is not printing, as the GOMC_FF_filename variable was not supplied a string for its name')

    else:
        print("******************************")
        print("")
        print('writing the GOMC force field file ')
        with open(GOMC_FF_filename, 'w') as data:
            if structure_1 != None:
                data.write("*  "+filename_0+' and '+ filename_1+
                           ' - created by mBuild using the on ' + str(date_time) +'\n') #
            else:
                data.write("*  " + filename_0 + ' - created by mBuild using the on ' + str(date_time) + '\n')
            if forcefield_files != None:
                data.write("*  " + 'parameters from the ' + str(forcefield_files) + ' force field(s) via MoSDef\n')
            elif forcefield_names != None:
                data.write("*  " + 'parameters from the '+str(forcefield_names)+' force field(s) via MoSDef\n')
            data.write("*  1-4 coulombic scaling = " + str(combined_1_4_Coul_dict_per_residue)+
                       ', and 1-4 LJ scaling = ' + str(combined_1_4_LJ_dict_per_residue)+'\n\n')
            data.write("*  "+'{:d} atoms\n'.format(len(structure_selection.atoms)))

            if atom_style in ['full', 'molecular']:
                data.write("*  "+'{:d} bonds\n'.format(len(bonds)))
                data.write("*  "+'{:d} angles\n'.format(len(angles)))
                data.write("*  "+'{:d} dihedrals\n'.format(len(dihedrals)))
                data.write("*  "+'{:d} impropers\n\n'.format(len(impropers)))

            data.write("*  "+'{:d} atom types\n'.format(len(set(types))))
            if atom_style in ['full', 'molecular']:
                if bonds:
                    data.write("*  "+'{:d} bond types\n'.format(len(set(unique_bond_types))))
                if angles:
                    data.write("*  "+'{:d} angle types\n'.format(len(set(unique_angle_types))))
                if dihedrals:
                    data.write("*  "+'{:d} dihedral types\n'.format(len(set(unique_dihedral_types))))
                if impropers:
                    data.write("*  "+'{:d} improper types\n'.format(len(set(unique_improper_types))))


            data.write('\n')


            data.write('\n*  Masses\n\n')
            data.write('! atom_types \tmass \t\t  atomTypeForceFieldName_ResidueName '
                       +'(i.e., atoms_type_per_utilized_FF)  \n')
            for atom_type,mass in mass_dict.items():
                mass_format = '*  {}\t\t{:.6f}\t! {}\n'
                data.write(mass_format.format(base10_to_base62_alph_num(atom_type),
                                              mass, unique_types[atom_type - 1]))


            # Bond coefficients
            if bonds:
                data.write('\n')
                data.write('BONDS * harmonic\n')
                data.write('!\n')
                data.write('!V(bond) = Kb(b - b0)**2\n')
                data.write('!\n')
                data.write('!Kb: kcal/mole/A**2\n')
                data.write('!b0: A\n')
                data.write('!Kb (kcal/mol) = Kb (K) * Boltz. const.; (9999999999 if no stretching)\n')
                data.write('!\n')



                if unit_style == 'real':
                    data.write('!atom_types \t Kb\tb0 \t\t  atoms_types_per_utilized_FF\n')
                elif unit_style == 'lj':
                    data.write('ERROR invalid option')
                for params,idx in unique_bond_types.items():
                    bond_format = '{}\t{}\t{}\t{}\t\t! {}\t{}\n'
                    if (fix_res_bonds_angles != None) and ((params[3] and  params[4]) in fix_res_bonds_angles ):
                        fix_bond_K_value = '999999999999'
                        data.write( bond_format.format(base10_to_base62_alph_num(atom_types_to_index_value_dict[params[2][0]+'_' + str(params[3])]),
                                                       base10_to_base62_alph_num(atom_types_to_index_value_dict[params[2][1]+'_' + str(params[4])]),
                                                       fix_bond_K_value, params[1],
                                                       params[2][0]+'_' + str(params[3]),
                                                       params[2][1]+'_' + str(params[4])))

                    else:
                        data.write( bond_format.format(base10_to_base62_alph_num(atom_types_to_index_value_dict[params[2][0] +'_' + str(params[3])]),
                                                       base10_to_base62_alph_num(atom_types_to_index_value_dict[params[2][1]+'_' + str(params[4])]),
                                                       params[0],params[1],
                                                       params[2][0]+'_' + str(params[3]),
                                                       params[2][1]+'_' + str(params[4])))


            # Angle coefficients
            if angles:
                if use_urey_bradleys:
                    data.write('\n!  Error: Urey Bradley terms detected but not written,'+
                               'since they are currently not compatible with GOMC\n')

                data.write('\nANGLES * harmonic\n')
                data.write('!\n')
                data.write('!V(angle) = Ktheta(Theta - Theta0)**2\n')
                data.write('!\n')
                data.write('!Ktheta: kcal/mole/rad**2\n')
                data.write('!Theta0: degrees\n')
                data.write('!\n')
                data.write('! Ktheta (kcal/mol) = Ktheta (K) * Boltz. const.\t\t\n')
                data.write('!\n')
                data.write('!atom_types \t\tKtheta\tTheta0\t\t\t  atoms_types_per_utilized_FF\n')
                for params,idx in unique_angle_types.items():

                    if (fix_res_bonds_angles != None) and ((params[4] and  params[5] and  params[6])
                                                           in fix_res_bonds_angles ):
                        fix_angle_K_value = '999999999999'
                        angle_format = '{}\t{}\t{}\t{}\t{:.5f}\t\t! {}\t{}\t{}\n'
                        data.write(angle_format.format(base10_to_base62_alph_num(atom_types_to_index_value_dict[params[3][0]+'_'+params[4]]),
                                                       base10_to_base62_alph_num(atom_types_to_index_value_dict[params[2]+'_'+params[5]]),
                                                       base10_to_base62_alph_num(atom_types_to_index_value_dict[params[3][1]+'_'+params[6]]),
                                                       fix_angle_K_value ,params[1],
                                                       params[3][0]+'_'+params[4],
                                                       params[2]+'_'+params[5],
                                                       params[3][1]+'_'+params[6]))

                    else:
                        data.write(
                            '{}\t{}\t{}\t{}\t{:.5f}\t\t! {}\t{}\t{}\n'.format(base10_to_base62_alph_num(atom_types_to_index_value_dict[params[3][0]+'_'+params[4]]),
                                                                              base10_to_base62_alph_num(atom_types_to_index_value_dict[params[2]+'_'+params[5]]),
                                                                              base10_to_base62_alph_num(atom_types_to_index_value_dict[params[3][1]+'_'+params[6]]),
                                                                              params[0], params[1],
                                                                              params[3][0]+'_'+params[4],
                                                                              params[2]+'_'+params[5],
                                                                              params[3][1]+'_'+params[6]))


            # Dihedral coefficients
            if dihedrals:
                if use_rb_torsions:
                    List_if_large_error_dihedral_overall = []

                    List_if_largest_error_abs_Values_for_dihedral_overall = []
                    List_dihedral_overall_error = []


                    List_if_abs_max_Values_for_dihedral_overall = []
                    List_Dihedral_atoms_all_dihedral_overall = []

                    data.write('\nDIHEDRALS * CHARMM\n')
                    data.write('!\n')
                    data.write('!V(dihedral) = Kchi(1 + cos(n(chi) - delta))\n')
                    data.write('!\n')
                    data.write('!Kchi: kcal/mole\n')
                    data.write('!n: multiplicity\n')
                    data.write('!delta: degrees\n')
                    data.write('!\n')
                    data.write('! Kchi (kcal/mol) = Kchi (K) * Boltz. const.\n')
                    data.write('! Boltzmann = 0.0019872041 kcal / (mol * K)\n')
                    data.write('!\n')
                    if unit_style == 'real':
                        data.write('!atom_types \t\t\tKchi\t\tn\tdelta\t\t  atoms_types_per_utilized_FF\n')

                    elif unit_style == 'lj':
                        data.write('ERROR invalid option')

                    for params,idx in unique_dihedral_types.items():
                        CHARMM_coeffs = RB_to_CHARMM(params[0],
                                                 params[1],
                                                 params[2],
                                                 params[3],
                                                 params[4],
                                                 params[5])

                        # check the error between the convertions of RB_tortions to CHARMM DIHEDRALS (start)
                        RB_to_CHARMM_Abs_diff = []
                        Pi = np.pi
                        Dihedral_steps =2*10**(-3)
                        Dihedral_range = 4*Pi
                        Dihedral_No_Steps = int(Dihedral_range/Dihedral_steps)+1


                        for i in range(0,  Dihedral_No_Steps+1):
                            t=i*Dihedral_steps

                            RB_dihedral_calc = params[0]  \
                                               + params[1] * (np.cos(t - Pi))** 1 \
                                               + params[2] * (np.cos(t - Pi)) ** 2 \
                                               + params[3] * (np.cos( t - Pi)) ** 3 \
                                               + params[4] * (np.cos( t - Pi)) ** 4 \
                                               + params[5] * (np.cos(t - Pi)) ** 5
                            """CHARMM_torsions 
                            = K0 * (1 + Cos[n0 * (t) - (d0)]) + K1 * (1 + Cos[n1 * (t) - (d1)]) + K2 * (
                                       # 1 + Cos[n2 * (t) - (d2)])
                            + K3 * (1 + Cos[n3 * (t) - (d3)]) + K4 * (1 + Cos[n4 * (t) - (d4)]) + K5 * (
                                       # 1 + Cos[n5 * (t) - (d5)])
                                       
                            = K0 + K1 * (1 + Cos[n1 * (t) - (d1)]) + K2 * (1 + Cos[n2 * (t) - (d2)])
                            + K3 * (1 + Cos[n3 * (t) - (d3)]) + K4 * (1 + Cos[n4 * (t) - (d4)]) + K5 * (
                                       # 1 + Cos[n5 * (t) - (d5)]). """

                            RB_to_CHARMM_calc =   CHARMM_coeffs[0, 0] * (1 + np.cos(int(CHARMM_coeffs[0, 1]) * (t)
                                                                                    - CHARMM_coeffs[0, 2]*Pi/180) )\
                                                  + CHARMM_coeffs[1, 0] * (1 + np.cos(int(CHARMM_coeffs[1, 1]) * (t)
                                                                                      - CHARMM_coeffs[1, 2]*Pi/180) ) \
                                                + CHARMM_coeffs[2, 0] * (1 + np.cos(int(CHARMM_coeffs[2, 1]) * (t)
                                                                                    - CHARMM_coeffs[2, 2]*Pi/180) ) \
                                                + CHARMM_coeffs[3, 0] * (1 + np.cos(int(CHARMM_coeffs[3, 1]) * (t)
                                                                                    - CHARMM_coeffs[3, 2]*Pi/180) ) \
                                                + CHARMM_coeffs[4, 0] * (1 + np.cos(int(CHARMM_coeffs[4, 1]) * (t)
                                                                                    - CHARMM_coeffs[4, 2]*Pi/180) ) \
                                                + CHARMM_coeffs[5, 0] * (1 + np.cos(int(CHARMM_coeffs[5, 1]) * (t)
                                                                                    - CHARMM_coeffs[5, 2]*Pi/180) )


                            RB_to_CHARMM_absolute_difference = np.absolute(RB_dihedral_calc-RB_to_CHARMM_calc)
                            RB_to_CHARMM_Abs_diff.append(RB_to_CHARMM_absolute_difference)


                        List_if_large_error_dihedral_iteration=[]
                        List_abs_max_dihedral_iteration = []



                        if max(RB_to_CHARMM_Abs_diff) > 10**(-10):
                            List_if_large_error_dihedral_iteration.append(1)
                            List_abs_max_dihedral_iteration.append(max(RB_to_CHARMM_Abs_diff))

                            List_if_large_error_dihedral_overall.append(1)
                            List_if_largest_error_abs_Values_for_dihedral_overall.append(max(RB_to_CHARMM_Abs_diff))
                            List_dihedral_overall_error.append(str(params[8])+', '+str(params[9])+', '
                                                                              +str(params[10])+', '+str(params[11]))



                        else:
                            List_if_large_error_dihedral_iteration.append(0)

                            List_if_abs_max_Values_for_dihedral_overall.append(max(RB_to_CHARMM_Abs_diff))
                            List_Dihedral_atoms_all_dihedral_overall.append(
                                str(params[8]) + ', ' + str(params[9]) + ', ' + str(params[10]) + ', ' + str(
                                    params[11]))

                        # **************************************
                        # check the error between the convertions of RB_tortions to CHARMM DIHEDRALS (end)
                        # **************************************
                        dihedral_format = '{}\t{}\t{}\t{}\t{:.6f}\t{}\t{}\t\t! {}\t{}\t{}\t{}\n'
                        data.write(dihedral_format.format(base10_to_base62_alph_num(atom_types_to_index_value_dict[ params[8]+ '_' + params[12]]),
                                                          base10_to_base62_alph_num(atom_types_to_index_value_dict[params[9]+ '_' + params[13]]),
                                                          base10_to_base62_alph_num(atom_types_to_index_value_dict[params[10]+ '_' + params[14]]),
                                                          base10_to_base62_alph_num(atom_types_to_index_value_dict[params[11]+ '_' + params[15]]),
                                                          CHARMM_coeffs[0,0],
                                                          int(CHARMM_coeffs[0,1]),
                                                          CHARMM_coeffs[0,2],
                                                          params[8]+ '_' + params[12],
                                                          params[9]+ '_' + params[13],
                                                          params[10]+ '_' + params[14],
                                                          params[11]+ '_' + params[15]))
                        data.write(dihedral_format.format(base10_to_base62_alph_num(atom_types_to_index_value_dict[params[8]+ '_' + params[12]]),
                                                          base10_to_base62_alph_num(atom_types_to_index_value_dict[params[9]+ '_' + params[13]]),
                                                          base10_to_base62_alph_num(atom_types_to_index_value_dict[params[10]+ '_' + params[14]]),
                                                          base10_to_base62_alph_num(atom_types_to_index_value_dict[params[11]+ '_' + params[15]]),
                                                          CHARMM_coeffs[1, 0],
                                                          int(CHARMM_coeffs[1, 1]),
                                                          CHARMM_coeffs[1, 2],
                                                          params[8]+ '_' + params[12],
                                                          params[9]+ '_' + params[13],
                                                          params[10]+ '_' + params[14],
                                                          params[11]+ '_' + params[15]))
                        data.write(dihedral_format.format(base10_to_base62_alph_num(atom_types_to_index_value_dict[params[8]+ '_' + params[12]]),
                                                          base10_to_base62_alph_num(atom_types_to_index_value_dict[params[9]+ '_' + params[13]]),
                                                          base10_to_base62_alph_num(atom_types_to_index_value_dict[params[10]+ '_' + params[14]]),
                                                          base10_to_base62_alph_num(atom_types_to_index_value_dict[params[11]+ '_' + params[15]]),
                                                          CHARMM_coeffs[2, 0],
                                                          int(CHARMM_coeffs[2, 1]),
                                                          CHARMM_coeffs[2, 2],
                                                          params[8]+ '_' + params[12],
                                                          params[9]+ '_' + params[13],
                                                          params[10]+ '_' + params[14],
                                                          params[11]+ '_' + params[15]))
                        data.write(dihedral_format.format(base10_to_base62_alph_num(atom_types_to_index_value_dict[params[8]+ '_' + params[12]]),
                                                          base10_to_base62_alph_num(atom_types_to_index_value_dict[params[9]+ '_' + params[13]]),
                                                          base10_to_base62_alph_num(atom_types_to_index_value_dict[params[10]+ '_' + params[14]]),
                                                          base10_to_base62_alph_num(atom_types_to_index_value_dict[params[11]+ '_' + params[15]]),
                                                          CHARMM_coeffs[3, 0],
                                                          int(CHARMM_coeffs[3, 1]),
                                                          CHARMM_coeffs[3, 2],
                                                          params[8]+ '_' + params[12],
                                                          params[9]+ '_' + params[13],
                                                          params[10]+ '_' + params[14],
                                                          params[11]+ '_' + params[15]))
                        data.write(dihedral_format.format(base10_to_base62_alph_num(atom_types_to_index_value_dict[params[8]+ '_' + params[12]]),
                                                          base10_to_base62_alph_num(atom_types_to_index_value_dict[params[9]+ '_' + params[13]]),
                                                          base10_to_base62_alph_num(atom_types_to_index_value_dict[params[10]+ '_' + params[14]]),
                                                          base10_to_base62_alph_num(atom_types_to_index_value_dict[params[11]+ '_' + params[15]]),
                                                          CHARMM_coeffs[4, 0],
                                                          int(CHARMM_coeffs[4, 1]),
                                                          CHARMM_coeffs[4, 2],
                                                          params[8]+ '_' + params[12],
                                                          params[9]+ '_' + params[13],
                                                          params[10]+ '_' + params[14],
                                                          params[11]+ '_' + params[15]))
                        data.write(dihedral_format.format(base10_to_base62_alph_num(atom_types_to_index_value_dict[params[8]+ '_' + params[12]]),
                                                          base10_to_base62_alph_num(atom_types_to_index_value_dict[params[9]+ '_' + params[13]]),
                                                          base10_to_base62_alph_num(atom_types_to_index_value_dict[params[10]+ '_' + params[14]]),
                                                          base10_to_base62_alph_num(atom_types_to_index_value_dict[params[11]+ '_' + params[15]]),
                                                          CHARMM_coeffs[5, 0],
                                                          int(CHARMM_coeffs[5, 1]),
                                                          CHARMM_coeffs[5, 2],
                                                          params[8]+ '_' + params[12],
                                                          params[9]+ '_' + params[13],
                                                          params[10]+ '_' + params[14],
                                                          params[11]+ '_' + params[15]))


                    if sum(List_if_large_error_dihedral_overall) > 0 :
                        List_if_largest_error_abs_Values_for_dihedral_overall_max = \
                            max(List_if_largest_error_abs_Values_for_dihedral_overall)
                        info_if_dihedral_error_too_large = '! WARNING: RB-torsion to CHARMM dihedral conversion error'\
                                                           +' is to large [error > 10^(-10)] \n' + \
                                                           '! WARNING: Maximum( '\
                                                           +'|(RB-torsion calc)-(CHARMM dihedral calc)| ) =  ' \
                                                           + str(List_if_largest_error_abs_Values_for_dihedral_overall_max)  \
                                                           +'\n'
                        data.write(info_if_dihedral_error_too_large)
                        print(info_if_dihedral_error_too_large)
                    else:
                        List_if_abs_max_Values_for_dihedral_overall_max = \
                            max(List_if_abs_max_Values_for_dihedral_overall)
                        info_if_dihedral_error_OK = '! RB-torsion to CHARMM dihedral conversion error is OK '\
                                                           +'[error <= 10^(-10)]\n' + \
                                                    '! Maximum( |(RB-torsion calc)-(CHARMM dihedral calc)| ) =  ' \
                                                    + str(List_if_abs_max_Values_for_dihedral_overall_max) + '\n'
                        data.write(info_if_dihedral_error_OK)
                        print(info_if_dihedral_error_OK)




                elif use_dihedrals:
                    data.write("Error not set up to use to use_dihedrals form for data input from the xml file")


            # Improper coefficients
            if impropers:

                data.write("Error GOMC is not currently able to use improper in its calculations")


            # Non-Bonded forces
            epsilons = np.array([atom.epsilon for atom in structure_selection.atoms]) / epsilon_conversion_factor
            sigmas = np.array([atom.sigma for atom in structure_selection.atoms]) / sigma_conversion_factor
            forcefields = [atom.type+'_'+atom.residue.name for atom in structure_selection.atoms]
            Residues = [atom.residue.name for atom in structure_selection.atoms]
            epsilon_dict = dict([(unique_types.index(atom_type)+1,epsilon)
                                 for atom_type,epsilon in zip(types,epsilons)])
            sigma_dict = dict([(unique_types.index(atom_type)+1,sigma) for atom_type,sigma in zip(types,sigmas)])
            LJ_1_4_dict = dict([(unique_types.index(atom_type) + 1, combined_1_4_LJ_dict_per_residue[Residues])
                                for atom_type, Residues in zip(types, Residues)])
            forcefield_dict = dict([(unique_types.index(atom_type) + 1, forcefield)
                                    for atom_type, forcefield in zip(types, forcefields)])

            # ensure all 1,4-coulombic scaling factors are the same
            coul_1_4_dict = dict(
                [(unique_types.index(atom_type) + 1, combined_1_4_Coul_dict_per_residue[Residues]) for
                 atom_type, Residues in zip(types, Residues)])
            coul_1_4_List =[]
            for p in coul_1_4_dict.values():
                coul_1_4_List.append(p)
            for c in range(0, len(coul_1_4_List)):
                if coul_1_4_List[c] != coul_1_4_List[0]:
                    return warn("ERROR: There are multiple 1,4-coulombic scaling factors.'+ "
                                "' GOMC will only accept a singular input for the 1,4-coulombic scaling factors")


            # Pair coefficients
            print('NBFIX_Mixing not used or no mixing used for the non-bonded potentials out')

            if non_bonded_type=='LJ':
                data.write('\n')
                data.write('NONBONDED\n')
                data.write('!\n')
                data.write('!V(Lennard-Jones) = Eps,i,j[(Rmin,i,j/ri,j)**12 - 2(Rmin,i,j/ri,j)**6]\n')
                data.write('!\n')

                if unit_style == 'real':
                    data.write('!atype \tignored\tepsilon \tRmin/2 \t\tignored\teps,1-4\t\tRmin/2,1-4\t\t'+
                               '  atom_type_per_utilized_FF\n')

                print('forcefield_dict = '+str(forcefield_dict))

                for idx, epsilon in epsilon_dict.items():
                    NB_format = '{}\t{:.2f}\t{:.9f}\t{:.11f}\t{:.2f}\t{:.9f}\t{:.11f}\t\t! {}\t{}\n'
                    data.write(NB_format.format(base10_to_base62_alph_num(idx), 0, -epsilon,
                                                sigma_dict[idx] * (2 ** (1 / 6)) / 2, 0,
                                                float(LJ_1_4_dict[idx])* (-epsilon),
                                                float(LJ_1_4_dict[idx])* sigma_dict[idx] * (2 ** (1 / 6)) / 2,
                                                forcefield_dict[idx],forcefield_dict[idx]))



            elif  non_bonded_type=='Mie':
                data.write("Error: Currenly the Mie potential is not supported in this MoSDeF GOMC parameter writer")
            else:
                data.write("Error: Currenly this potential is not supported in this MoSDeF GOMC parameter writer")


            data.write('\n!\n')

            if forcefield_files != None:
                data.write('!NBFIX not used for the '+str(forcefield_files)+' file type(s)')
            elif forcefield_names != None:
                data.write('!NBFIX not used for the '+str(forcefield_names)+' file type(s)')



    # **********************************
    #**********************************
    # FF writer (end)
    # **********************************
    # **********************************








    # **********************************
    #**********************************
    # psf writer (start)
    # **********************************
    # **********************************


    print("******************************")
    print("")
    print('write_charmm_psf file is running')
    vmd = True

    date_time = datetime.datetime.today()

    if forcefield_files != None:
        print('write_charmm_psf: forcefield_files = ' + str(forcefield_files) + ', ' + 'residues = ' + str(residues))
    elif forcefield_names != None:
        print('write_charmm_psf: forcefield_names = ' + str(forcefield_names) + ', ' + 'residues = ' + str(residues))


    print("******************************")
    print("")


    if structure_1 !=None:
        list_of_structures = [structure_0_FF, structure_1_FF]
        list_of_file_names = [filename_0, filename_1]
        stuct_only = [structure_0_FF, structure_1_FF]
    else:
        list_of_structures = [structure_0_FF]
        list_of_file_names = [filename_0]
        stuct_only = [structure_0_FF]

    for q in range(0, len(list_of_structures)):
        stuct_iteration = list_of_structures[q]
        file_name_iteration = list_of_file_names[q]
        dest = str(file_name_iteration)+'.psf'
        stuct_only_iteration =stuct_only[q]
        # Lammps syntax depends on the functional form
        # Infer functional form based on the properties of the stuct_iteration
        if detect_forcefield_style:
            # Check angles
            if len(stuct_iteration.urey_bradleys) > 0:
                print("Urey bradley terms detected")
                data.write("Urey bradley terms detected")
                data.write("Warning: GOMC does no support the Urey-Bradley terms")
                use_urey_bradleys = True
            else:
                print("No urey bradley terms detected")
                use_urey_bradleys = False

            # Check dihedrals
            if len(stuct_iteration.rb_torsions) > 0:
                print("RB Torsions detected, will converted to CHARMM Dihedrals")
                use_rb_torsions = True
                dihedral_rb_torsions_stuct_iteration_list = stuct_iteration.rb_torsions
            else:
                use_rb_torsions = False

            if len(stuct_iteration.dihedrals) > 0:
                print("Charmm dihedrals detected, so CHARMM Dihedrals will remain")
                use_dihedrals = True
                dihedral_rb_torsions_stuct_iteration_list = stuct_iteration.dihedrals
            else:
                use_dihedrals = False
        if (use_rb_torsions == False) and (use_dihedrals == False):
            dihedral_rb_torsions_stuct_iteration_list = []
        if use_rb_torsions and use_dihedrals:
            raise ValueError("Multiple dihedral styles detected, check your "
                             "Forcefield XML and stuct_iteration")

        # Check impropers
        for dihedral in stuct_iteration.dihedrals:
            if dihedral.improper:
                raise ValueError("Error: Amber-style impropers are currently not supported")

        dihedral_impropers_stuct_iteration_list = stuct_iteration.impropers


        # psf printing (start)


        own_handle = False

        residue_data_list = []
        residue_names_list = []
        for k, atom in enumerate(stuct_only_iteration.atoms):
            residue_data_list.append(str(atom.residue))
            residue_names_list.append(atom.residue.name)

        unique_residue_data_dict = {}
        unique_residue_data_list = []
        residue_data_name_list = []

        for m, residue in enumerate(stuct_only_iteration.residues):
            unique_residue_data_list.append(str(stuct_only_iteration.residues[m]))
            unique_residue_data_dict.update({unique_residue_data_list[m]: m + 1})
            residue_data_name_list.append(stuct_only_iteration.residues[m].name)


        Max_Residue_No = 9999
        No_last_values_res_name = 3

        Res_No_iteration_corrected_List = []
        residue_ID_list = []
        for f, atom in enumerate(stuct_only_iteration.atoms):

            residue_ID_int = int(unique_residue_data_dict[residue_data_list[f]])
            Res_ID_adder = int((residue_ID_int % Max_Residue_No) % (Max_Residue_No))
            if int(Res_ID_adder) == 0:
                Res_No_iteration_corrected = int(Max_Residue_No)
            else:
                Res_No_iteration_corrected = Res_ID_adder

            Res_No_iteration_corrected_List.append(Res_No_iteration_corrected)
            residue_ID_list.append(residue_ID_int)

        # Index the atoms and residues TODO delete
        if isinstance(dest, string_types):
            own_handle = True
            dest = genopen(dest, 'w')

        atmfmt1 = ('%8s %-4s %-4s %-4s %-4s %4s %10.6f %13.4f' + 11 * ' ')
        intfmt = '%8s'  # For pointers

        dest.write('PSF ')  # BC Changed from 'PSF CHEQ ' to 'PSF'
        dest.write('\n\n')
        if isinstance(stuct_iteration.title, string_types):
            dest.write(intfmt % 1 + ' !NTITLE\n')
            dest.write(' REMARKS this file ' + file_name_iteration + ' - created by mBuild/foyer using the' + '\n')
            if forcefield_files != None:
                dest.write(' REMARKS parameters from the ' + str(forcefield_files) + ' force field via MoSDef\n')
            elif forcefield_names != None:
                dest.write(' REMARKS parameters from the ' + str(forcefield_names) + ' force field via MoSDef\n')
            dest.write(' REMARKS created on ' + str(date_time) + '\n')
            dest.write('%s\n\n' % stuct_iteration.title)
        else:
            dest.write(intfmt % len(stuct_iteration.title) + ' !NTITLE\n')
            dest.write(' REMARKS this file ' + file_name_iteration + ' - created by mBuild/foyer using the' + '\n')
            if forcefield_files != None:
                dest.write(' REMARKS parameters from the ' + str(forcefield_files) + ' force field via MoSDef\n')
            elif forcefield_names != None:
                dest.write(' REMARKS parameters from the ' + str(forcefield_names) + ' force field via MoSDef\n')
            dest.write(' REMARKS created on ' + str(date_time) + '\n')
            dest.write('\n'.join(stuct_iteration.title) + '\n\n')


        # This converts the atom name in the GOMC psf and pdb files to unique atom names
        unique_Individual_atom_names_dict, \
        Individual_atom_names_List, \
        Missing_Bead_to_atom_name =unique_atom_naming(stuct_only_iteration , residue_ID_list, residue_names_list,
                                                      Bead_to_atom_name_dict=Bead_to_atom_name_dict)

        # Now time for the atoms
        dest.write(intfmt % len(stuct_iteration.atoms) + ' !NATOM\n')
        # atmfmt1 is for CHARMM format (i.e., atom types are integers)
        for i, atom in enumerate(stuct_iteration.atoms):
            typ = base10_to_base62_alph_num(atom_types_to_index_value_dict[atom.type+'_'+atom.residue.name])
            fmt = atmfmt1
            segid = atom.residue.segid or 'SYS'


            atmstr = fmt % (i + 1, segid,
                            Res_No_iteration_corrected_List[i],
                            str(residue_names_list[i])[:No_last_values_res_name], Individual_atom_names_List[i], typ,
                            atom.charge, atom.mass)

            if hasattr(atom, 'props'):
                dest.write(atmstr + '   '.join(atom.props) + '\n')
            else:
                dest.write('%s\n' % atmstr)
        dest.write('\n')

        # Bonds
        dest.write(intfmt % len(stuct_iteration.bonds) + ' !NBOND: bonds\n')
        for i, bond in enumerate(stuct_iteration.bonds):

            dest.write((intfmt * 2) % (bond.atom1.idx + 1, bond.atom2.idx + 1))
            if i % 4 == 3:  # Write 4 bonds per line
                dest.write('\n')
        # See if we need to terminate
        if len(stuct_iteration.bonds) % 4 != 0 or len(stuct_iteration.bonds) == 0:
            dest.write('\n')
        dest.write('\n')

        # Angles
        dest.write(intfmt % len(stuct_iteration.angles) + ' !NTHETA: angles\n')
        for i, angle in enumerate(stuct_iteration.angles):

            dest.write((intfmt * 3) % (angle.atom1.idx + 1, angle.atom2.idx + 1,
                                       angle.atom3.idx + 1)
                       )
            if i % 3 == 2:  # Write 3 angles per line
                dest.write('\n')
        # See if we need to terminate
        if len(stuct_iteration.angles) % 3 != 0 or len(stuct_iteration.angles) == 0:
            dest.write('\n')
        dest.write('\n')

        # Dihedrals
        # impropers need to be split off in the "improper" section.
        # PSF files need to have each dihedral listed *only* once. So count the
        # number of unique dihedrals
        nnormal = 0
        torsions = set()
        # for dih in stuct_iteration.dihedrals:
        for dih in dihedral_rb_torsions_stuct_iteration_list:
            if dih.improper: continue
            a1, a2, a3, a4 = dih.atom1, dih.atom2, dih.atom3, dih.atom4
            if (a1, a2, a3, a4) in torsions or (a4, a3, a2, a1) in torsions:
                continue
            nnormal += 1
            torsions.add((a1, a2, a3, a4))
        nimprop = sum(1 for dih in stuct_iteration.dihedrals if dih.improper)
        dest.write(intfmt % nnormal + ' !NPHI: dihedrals\n')
        torsions = set()
        c = 0
        # for dih in stuct_iteration.dihedrals:
        for dih in dihedral_rb_torsions_stuct_iteration_list:
            if dih.improper: continue
            a1, a2, a3, a4 = dih.atom1, dih.atom2, dih.atom3, dih.atom4
            if (a1, a2, a3, a4) in torsions or (a4, a3, a2, a1) in torsions:
                continue

            dest.write((intfmt * 4) % (a1.idx + 1, a2.idx + 1, a3.idx + 1, a4.idx + 1))
            torsions.add((a1, a2, a3, a4))
            if c % 2 == 1:  # Write 2 dihedrals per line
                dest.write('\n')
            c += 1
        # See if we need to terminate
        if nnormal % 2 != 0 or nnormal == 0:
            dest.write('\n')
        dest.write('\n')
        # Impropers
        nimprop += len(dihedral_impropers_stuct_iteration_list)
        dest.write(intfmt % (nimprop) + ' !NIMPHI: impropers\n')

        def improp_gen(stuct_iteration):
            # for imp in stuct_iteration.impropers:
            for imp in dihedral_impropers_stuct_iteration_list:
                yield (imp.atom1, imp.atom2, imp.atom3, imp.atom4)
            # for dih in stuct_iteration.dihedrals:
            for dih in dihedral_rb_torsions_stuct_iteration_list:
                if dih.improper:
                    yield (dih.atom1, dih.atom2, dih.atom3, dih.atom4)

        for i, (a1, a2, a3, a4) in enumerate(improp_gen(stuct_iteration)):

            dest.write((intfmt * 4) % (a1.idx + 1, a2.idx + 1, a3.idx + 1, a4.idx + 1))
            if i % 2 == 1:  # Write 2 dihedrals per line
                dest.write('\n')
        # See if we need to terminate
        if nimprop % 2 != 0 or nimprop == 0:
            dest.write('\n')
        dest.write('\n')
        # Donor section
        dest.write(intfmt % len(stuct_iteration.donors) + ' !NDON: donors\n')
        for i, don in enumerate(stuct_iteration.donors):

            dest.write((intfmt * 2) % (don.atom1.idx + 1, don.atom2.idx + 1))
            if i % 4 == 3:  # 4 donors per line
                dest.write('\n')
        if len(stuct_iteration.donors) % 4 != 0 or len(stuct_iteration.donors) == 0:
            dest.write('\n')
        dest.write('\n')
        # Acceptor section
        dest.write(intfmt % len(stuct_iteration.acceptors) + ' !NACC: acceptors\n')
        for i, acc in enumerate(stuct_iteration.acceptors):

            dest.write((intfmt * 2) % (acc.atom1.idx + 1, acc.atom2.idx + 1))
            if i % 4 == 3:  # 4 donors per line
                dest.write('\n')
        if len(stuct_iteration.acceptors) % 4 != 0 or len(stuct_iteration.acceptors) == 0:
            dest.write('\n')
        dest.write('\n')
        # NNB section ??
        dest.write(intfmt % 0 + ' !NNB\n\n')
        for i in range(len(stuct_iteration.atoms)):
            dest.write(intfmt % 0)
            if i % 8 == 7:  # Write 8 0's per line
                dest.write('\n')
        if len(stuct_iteration.atoms) % 8 != 0: dest.write('\n')
        dest.write('\n')
        # Group section
        try:
            nst2 = stuct_iteration.groups.nst2
        except AttributeError:
            nst2 = 0
        dest.write((intfmt * 2) % (len(stuct_iteration.groups) or 1, nst2))
        dest.write(' !NGRP \n')  # Changed from ' !NGRP NST2\n' to ' !NGRP \n'
        if stuct_iteration.groups:
            for v, gp in enumerate(stuct_iteration.groups):

                dest.write((intfmt * 3) % (gp.atom.idx, gp.type, gp.move))
                if v % 3 == 2: dest.write('\n')
            if len(stuct_iteration.groups) % 3 != 0 or len(stuct_iteration.groups) == 0:
                dest.write('\n')
        else:
            typ = 1 if abs(sum(a.charge for a in stuct_iteration.atoms)) < 1e-4 else 2
            dest.write((intfmt * 3) % (0, typ, 0))
            dest.write('\n')
        dest.write('\n')

        # The next two sections are never found in VMD prmtops...
        if not vmd:
            # Molecule section; first set molecularity
            set_molecules(stuct_iteration.atoms)
            mollist = [a.marked for a in stuct_iteration.atoms]
            dest.write(intfmt % max(mollist) + ' !MOLNT\n')
            for v, atom in enumerate(stuct_iteration.atoms):
                dest.write(intfmt % atom.marked)
                if v % 8 == 7: dest.write('\n')
            if len(stuct_iteration.atoms) % 8 != 0: dest.write('\n')
            dest.write('\n')
            # NUMLP/NUMLPH section
            dest.write((intfmt * 2) % (0, 0) + ' !NUMLP NUMLPH\n')
            dest.write('\n')

            # CMAP section

            dest.write(intfmt % len(stuct_iteration.cmaps) + ' !NCRTERM: cross-terms\n')
            for v, cmap in enumerate(stuct_iteration.cmaps):
                dest.write((intfmt * 8) % (cmap.atom1.idx + 1, cmap.atom2.idx + 1,
                                           cmap.atom3.idx + 1, cmap.atom4.idx + 1,
                                           cmap.atom2.idx + 1, cmap.atom3.idx + 1,
                                           cmap.atom4.idx + 1, cmap.atom5.idx + 1)
                           )
                dest.write('\n')

        dest.write('\n')
        # Done!
        # If we opened our own handle, close it
        if own_handle:
            dest.close()
    # **********************************
    # **********************************
    # psf writer (end)
    # **********************************
    # **********************************







    # **********************************
    # **********************************
    # pdb writer (start)
    # **********************************
    # **********************************

    print("******************************")
    print("")
    print('write_charmm_pdb file is running')
    print('write_charmm_pdb: residues == ' + str(residues))
    print('fix_residue = ' + str(fix_residue))
    print('fix_residue_in_box = ' + str(fix_residue_in_box))
    print('Bead_to_atom_name_dict = ' + str(Bead_to_atom_name_dict))



    if fix_residue is None and fix_residue_in_box is None:
        print('INFORMATION: No atoms are fixed in this pdb file for the GOMC simulation engine. ')
    else:
        warn('Some atoms are fixed in this pdb file for the GOMC simulation engine. ')

    print("******************************")
    print("")



    if structure_1 !=None:
        list_of_structures = [structure_0_FF, structure_1_FF]
        list_of_file_names = [filename_0, filename_1]
        stuct_only = [structure_0_FF, structure_1_FF]
    else:
        list_of_structures = [structure_0_FF]
        list_of_file_names = [filename_0]
        stuct_only = [structure_0_FF]

    for q in range(0, len(list_of_structures)):
        file_name_iteration = list_of_file_names[q]
        dest = str(file_name_iteration)+'.pdb'
        stuct_only_iteration =stuct_only[q]

        dest = genopen(dest, 'w')


        unique_residue_data_dict = {}
        unique_residue_data_list = []
        residue_data_name_list = []

        for m, residue in enumerate( stuct_only_iteration.residues):
            unique_residue_data_list.append(str( stuct_only_iteration.residues[m]))
            unique_residue_data_dict.update({ unique_residue_data_list[m]: m+1})
            residue_data_name_list.append( stuct_only_iteration.residues[m].name)

        for n in range(0, len(residue_data_name_list)):
            if residue_data_name_list[n] not in  residues:
                return warn('Error: Please specifiy all residues (residues) in a list')


        residue_data_list = []
        for k, atom in enumerate( stuct_only_iteration.atoms):
            residue_data_list.append(str(atom.residue))

        if (fix_residue != None) and (fix_residue_in_box != None):
            for n in range(0,len(fix_residue)):
                if fix_residue[n] in fix_residue_in_box:
                    return warn("ERROR: residue type can not be specified to both fix_residue and fix_residue_in_box")

        residue_names_list =[]
        fix_atoms_list = []
        for k, atom in enumerate( stuct_only_iteration.atoms):
            residue_names_list.append(atom.residue.name)
            if (fix_residue != None) and (atom.residue.name  in fix_residue) :
                    beta_iteration = 1.00
            elif (fix_residue_in_box != None) and (atom.residue.name  in fix_residue_in_box):
                beta_iteration = 2.00
            else:
                beta_iteration = 0.00
            fix_atoms_list.append(beta_iteration)


        if  stuct_only_iteration.box is not None:
            dest.write('CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f %-11s%4s\n' % (
                 stuct_only_iteration.box[0],  stuct_only_iteration.box[1],
                 stuct_only_iteration.box[2],  stuct_only_iteration.box[3],
                 stuct_only_iteration.box[4],  stuct_only_iteration.box[5],
                 stuct_only_iteration.space_group, ''))
        if  stuct_only_iteration.symmetry is not None:
            fmt = '%d%4d%10.6f%10.6f%10.6f%15.5f\n'
            for index, arr in enumerate( stuct_only_iteration.symmetry.data):
                arr_list = [1 + index % 3, 1 + index // 3] + arr.tolist()
                symm_line = "REMARK 290   SMTRY" + fmt % tuple(arr_list)
                dest.write(symm_line)

        coords = stuct_only_iteration.get_coordinates('all')
        if coords is None:
            raise ValueError('Cannot write PDB file with no coordinates')
        # Create a function to process each atom and return which one we want
        # to print, based on our alternate location choice


        atomrec = ('ATOM  %5s %-4s%1s%-4s%1s%4d%1s   %8.3f%8.3f%8.3f%6.2f'
                       '%6.2f      %-4s%2s%-2s\n')

        pa_altloc_List = []
        res_chain_List = []
        res_insertion_code_List = []
        x_List = []
        y_List = []
        z_List = []
        pa_occupancy_List = []
        pa_bfactor_List = []
        Element_List = []

        # lock occupany factor at 1 (instead of: atom.occupancy)
        locked_occupany_factor = 1.00

        Max_No_atoms_in_base10 = 99999  # 99,999 for atoms in psf/pdb

        Max_Residue_No = 9999
        No_last_values_res_name = 3

        Res_No_iteration_corrected_List =[]
        residue_ID_list = []
        for i, atom in enumerate( stuct_only_iteration.atoms):

            residue_ID_int = int(unique_residue_data_dict[residue_data_list[i]])
            Res_ID_adder = int((residue_ID_int % Max_Residue_No) % (Max_Residue_No))
            if int(Res_ID_adder) == 0:
                Res_No_iteration_corrected_List.append(int(Max_Residue_No))
            else:
                Res_No_iteration_corrected_List.append(Res_ID_adder)

            residue_ID_list.append(residue_ID_int)

        # This converts the atom name in the GOMC psf and pdb files to unique atom names
        unique_Individual_atom_names_dict, \
        Individual_atom_names_List, \
        Missing_Bead_to_atom_name = unique_atom_naming(stuct_only_iteration, residue_ID_list, residue_names_list,
                                                       Bead_to_atom_name_dict=Bead_to_atom_name_dict)

        for  model, coord in enumerate(coords):


            if coords.shape[0] > 1:
                dest.write('MODEL      %5d\n' % (model+1))
            for res in  stuct_only_iteration.residues:
                atoms = sorted(res.atoms, key=lambda atom: atom.number)
                segid = ''
                for atom in atoms:
                    pa, others, (x, y, z) = print_atoms(atom, coord)
                    # Figure out the serial numbers we want to print


                    pa_altloc_List.append(pa.altloc)
                    res_chain_List.append(res.chain[-1:])
                    res_insertion_code_List.append(res.insertion_code[:1])
                    x_List.append(x)
                    y_List.append(y)
                    z_List.append(z)
                    pa_occupancy_List.append(pa.occupancy)
                    pa_bfactor_List.append(pa.bfactor)
                    Element_List.append(Element[pa.atomic_number].upper())


            for v, atom in enumerate( stuct_only_iteration.atoms):

                if v + 1 > Max_No_atoms_in_base10:
                    adj_atom_number = base10_to_base16_alph_num(v + 1)
                else:
                    adj_atom_number = v + 1

                dest.write(atomrec % (adj_atom_number, Individual_atom_names_List[v], pa_altloc_List[v],
                                      str(residue_names_list[v])[:No_last_values_res_name], res_chain_List[v],
                                      Res_No_iteration_corrected_List[v],
                                      res_insertion_code_List[v], x_List[v], y_List[v], z_List[v],
                                      locked_occupany_factor, fix_atoms_list[v], segid, Element_List[v], ''))


            dest.write("%-80s\n" % "END")

            if coords.shape[0] > 1:
                dest.write('ENDMDL\n')


        own_handle = True
        if own_handle:
            dest.close()






        # **********************************
        # **********************************
        # pdb writer (end)
        # **********************************
        # **********************************
