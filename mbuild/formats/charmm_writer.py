import os
import datetime
import numpy as np

from collections import OrderedDict
from warnings import warn
from parmed.utils.io import genopen

from parmed.periodic_table import Element
from parmed.utils.six.moves import range

from mbuild.compound import Compound
from mbuild import Box
from mbuild.utils.conversion import RB_to_CHARMM
from mbuild.utils.sorting import natural_sort
from mbuild.utils.conversion import base10_to_base16_alph_num
from mbuild.utils.conversion import base10_to_base26_alph
from mbuild.utils.conversion import base10_to_base52_alph
from mbuild.utils.conversion import base10_to_base62_alph_num
from mbuild.utils.specific_ff_to_residue import specific_ff_to_residue


def generate_impropers_for_PSF(stucture,
                               Dihedrals):

    for improper_iteration in stucture.impropers:
        yield (improper_iteration.atom1, improper_iteration.atom2,
               improper_iteration.atom3, improper_iteration.atom4)

    for dihedral_iteration in Dihedrals:
        if dihedral_iteration.improper:
            yield (dihedral_iteration.atom1, dihedral_iteration.atom2,
                   dihedral_iteration.atom3, dihedral_iteration.atom4)

def _get_bond_type_key(bond,
                       sigma_conversion_factor,
                       epsilon_conversion_factor):
    """Get the bond_type key for a bond"""
    return (
        round(bond.type.k * (sigma_conversion_factor ** 2 / epsilon_conversion_factor), 8),
        round(bond.type.req / sigma_conversion_factor, 8),
        tuple(sorted((bond.atom1.type, bond.atom2.type))),
        bond.atom1.residue.name, bond.atom2.residue.name
     )


def _get_angle_type_key(angle,
                        sigma_conversion_factor,
                        epsilon_conversion_factor):
    """Get the angle_type key for an angle"""
    return (
        round(angle.type.k*(sigma_conversion_factor**2/epsilon_conversion_factor), 8),
        round(angle.type.theteq, 8),
        angle.atom2.type,
        tuple(sorted((angle.atom1.type, angle.atom3.type))),
        angle.atom1.residue.name, angle.atom2.residue.name,
        angle.atom3.residue.name
    )


def _get_dihedral_rb_torsion_key(dihedral,
                                 epsilon_conversion_factor):
    lj_unit = 1 / epsilon_conversion_factor
    return (
        round(dihedral.type.c0*lj_unit, 8),
        round(dihedral.type.c1*lj_unit, 8),
        round(dihedral.type.c2*lj_unit, 8),
        round(dihedral.type.c3*lj_unit, 8),
        round(dihedral.type.c4*lj_unit, 8),
        round(dihedral.type.c5*lj_unit, 8),
        round(dihedral.type.scee, 1),
        round(dihedral.type.scnb, 1),
        dihedral.atom1.type,
        dihedral.atom2.type,
        dihedral.atom3.type,
        dihedral.atom4.type,
        dihedral.atom1.residue.name,
        dihedral.atom2.residue.name,
        dihedral.atom3.residue.name,
        dihedral.atom4.residue.name
    )


def _get_improper_type_key(improper,
                           epsilon_conversion_factor):
    lj_unit = 1 / epsilon_conversion_factor
    return (
        round(improper.type.psi_k * lj_unit, 8),
        round(improper.type.psi_eq, 8),
        improper.atom1.type,
        improper.atom2.type,
        improper.atom3.type,
        improper.atom4.type,
        improper.atom1.residue.name,
        improper.atom2.residue.name,
        improper.atom3.residue.name,
        improper.atom4.residue.name
    )


def _get_unique_bond_types(structure,
                           sigma_conversion_factor,
                           epsilon_conversion_factor):
    unique_bond_set = set()
    for bond in structure.bonds:
        unique_bond_set.add(
            _get_bond_type_key(bond, sigma_conversion_factor, epsilon_conversion_factor)
        )

    return {bond_key: i+1 for i, bond_key in enumerate(unique_bond_set)}


def _get_unique_angle_types(structure,
                            sigma_conversion_factor,
                            epsilon_conversion_factor):
    unique_angle_set = set()
    for angle in structure.angles:
        unique_angle_set.add(
            _get_angle_type_key(angle, sigma_conversion_factor, epsilon_conversion_factor)
        )

    return {angle_key: i+1 for i, angle_key in enumerate(unique_angle_set)}


def _get_unique_rb_torsion_types(structure,
                                 epsilon_conversion_factor):
    unique_dihedral_set = set()
    for dihedral in structure.rb_torsions:
        unique_dihedral_set.add(
            _get_dihedral_rb_torsion_key(dihedral, epsilon_conversion_factor)
        )

    return {torsion_key: i + 1 for i, torsion_key in enumerate(unique_dihedral_set)}


def _get_unique_improper_types(structure,
                               epsilon_conversion_factor):
    unique_improper_set = set()
    for improper in structure.impropers:
        unique_improper_set.add(_get_improper_type_key(improper, epsilon_conversion_factor))

    return {improper_key: i + 1 for i, improper_key in enumerate(unique_improper_set)}


def _get_bond_types(structure,
                    sigma_conversion_factor,
                    epsilon_conversion_factor):
    """Get a list of unique bond_types given a parmed structure

    Parameters
    ----------
    structure: parmed.Structure
        The parmed structure
    #ToDO: What is this intended for?
    #FIXME: Improve this docstring
    sigma_conversion_factor: float
        The sigma conversion factor
    epsilon_conversion_factor: float
        The epsilon conversion factor

    Returns
    -------

    """
    unique_bond_types = _get_unique_bond_types(structure, sigma_conversion_factor, epsilon_conversion_factor )

    bond_types = [
        unique_bond_types[_get_bond_type_key(bond, sigma_conversion_factor, epsilon_conversion_factor)]
        for bond in structure.bonds
    ]

    unique_bond_check_dict = {}
    for i_value_bond, i_key_bond in unique_bond_types.items():
        i_value_duplicated = False
        for j_value_bond, j_key_bond in unique_bond_types.items():
            j_value_bond_reorder = (j_value_bond[0], j_value_bond[1],
                                    j_value_bond[2][0], j_value_bond[2][0],
                                    j_value_bond[3], j_value_bond[4])

            if i_value_bond == j_value_bond_reorder:
                i_value_duplicated = True
                if i_value_bond[2][0] > j_value_bond[2][0]:
                    unique_bond_check_dict.update({j_value_bond: len(unique_bond_check_dict)})
                else:
                    unique_bond_check_dict.update({i_value_bond: len(unique_bond_check_dict)})

            if i_value_duplicated == False:
                unique_bond_check_dict.update({i_value_bond: len(unique_bond_check_dict)})

    unique_bond_types = OrderedDict([(y, x) for y, x in unique_bond_check_dict.items()])

    return bond_types, unique_bond_types


def _get_angle_types(structure,
                     sigma_conversion_factor,
                     epsilon_conversion_factor,
                     use_urey_bradleys=False):
    if use_urey_bradleys:
        warn('CRITICAL WARNING:  Urey-Bradleys are not available in the current '
             'version of this psf, pdb, and GOMC writer.')
        return None, None
    else:
        unique_angle_types = _get_unique_angle_types(structure, sigma_conversion_factor, epsilon_conversion_factor )

        angle_types = [unique_angle_types[_get_angle_type_key( angle, sigma_conversion_factor,
                                                               epsilon_conversion_factor
                                                               )] for angle in structure.angles]

    unique_angle_check_dict = {}
    for i_value_ang, i_key_ang in unique_angle_types.items():
        i_value_duplicated = False
        for j_value_ang, j_key_ang in unique_angle_types.items():
            j_value_ang_reorder = (j_value_ang[0], j_value_ang[1],
                                   j_value_ang[2], j_value_ang[3][0], j_value_ang[3][1],
                                     j_value_ang[4], j_value_ang[5], j_value_ang[6])

            if i_value_ang == j_value_ang_reorder:
                i_value_duplicated = True
                if i_value_ang[2] > j_value_ang[2]:
                    unique_angle_check_dict.update({j_value_ang: len(unique_angle_check_dict) })
                else:
                    unique_angle_check_dict.update({i_value_ang: len(unique_angle_check_dict) })

        if not i_value_duplicated:
            unique_angle_check_dict.update({i_value_ang: len(unique_angle_check_dict)})

    unique_angle_types = OrderedDict([(y, x) for y, x in unique_angle_check_dict.items()])

    return angle_types, unique_angle_types


def _get_dihedral_types(structure,
                        use_rb_torsions,
                        use_dihedrals,
                        epsilon_conversion_factor):
    if use_rb_torsions:
        unique_dihedral_types = _get_unique_rb_torsion_types(structure, epsilon_conversion_factor)

        dihedral_types = [
            unique_dihedral_types[
                _get_dihedral_rb_torsion_key(dihedral, epsilon_conversion_factor)]
            for dihedral in structure.rb_torsions
        ]

    elif use_dihedrals:
        warn('CRITICAL WARNING: Using the charmm style and impropers is not '
             'available in the current version of this psf, pdb, and GOMC writer.')
        return None, None

    unique_dihedral_check_dict = OrderedDict()
    for i_value_dihed, i_key_dihed in unique_dihedral_types.items():
        i_value_duplicated = False
        for j_value_dihed, j_key_dihed in unique_dihedral_types.items():
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
    unique_improper_types = _get_unique_improper_types(structure, epsilon_conversion_factor)

    improper_types = [
        unique_improper_types[
            _get_improper_type_key(improper, epsilon_conversion_factor)]
                for improper in structure.impropers
    ]

    # If impropers are added to GOMC, add the atom sorter for the unique combinations here

    return improper_types, unique_improper_types


def unique_atom_naming(structure, residue_ID_list, residue_names_list, bead_to_atom_name_dict=None):

    """ Outputs
        unique_Individual_atom_names_dict : dictionary
            All the unique atom names compiled into a dictionary.
        Individual_atom_names_List : list, in sequential  order
            The atom names for every atom in the system
        Missing_Bead_to_atom_name  : list, in sequential  order
            The bead names of any atoms beads that did not have a name specificed to them
            via the bead_to_atom_name_dict

    Parameters
    ----------
    structure : compound object
    residue_ID_list : list, in sequential  order
            The residue ID for every atom in the system
    residue_names_list  : list, in sequential  order
        The atom names for every atom in the system
    bead_to_atom_name_dict: dictionary ; optional, default =None
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
                if bead_to_atom_name_dict != None and (str(atom.name) in bead_to_atom_name_dict) == True:
                    if len(bead_to_atom_name_dict[str(atom.name)]) > 2:
                        text_to_write = ('ERROR: only enter atom names that have 2 or less digits' +
                                         ' in the Bead to atom naming dictionary (bead_to_atom_name_dict).')
                        warn(text_to_write)
                        return None, None, None
                    else:
                        atom_name_value = bead_to_atom_name_dict[str(atom.name)]
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
                    warn(text_to_write)
                    return None, None, None
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
        warn("NOTE: All bead names were not found in the Bead to atom naming dictionary (bead_to_atom_name_dict) ")

    return unique_Individual_atom_names_dict, Individual_atom_names_List, Missing_Bead_to_atom_name




# Currently the NBFIX is disabled as since only the OPLS and TRAPPE force fields are currently supported
class Charmm:
    def __init__(self, structure_box_0, filename_box_0, structure_box_1 = None, filename_box_1= None,
                 non_bonded_type='LJ', forcefield_selection = None,  residues=None,
                 detect_forcefield_style=True, fix_res_bonds_angles = None, bead_to_atom_name_dict=None,
                 fix_residue=None, fix_residue_in_box=None,  FF_filename= None,
                 reorder_res_in_pdb_psf =False, box_0 = None, box_1 = None  , **kwargs):

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
        structure_box_0 : compound object
        filename_box_0 : str
            Path of the output file for structure_box_0
        structure_box_1 : compound object number 2, optional
            (Ex: for GCMC or GEMC simulations which have mulitiple simulation boxes)
        filename_box_1 : str , optional
            Path of the output file for structure_box_1 (Ex: for GCMC or GEMC simulations
            which have mulitiple simulation boxes)
        non_bonded_type : str. optional, default = 'LJ' (i.e., Lennard-Jones )
                Specify the type of non-bonded potential for the GOMC force field files.
                Note: Currently, on the 'LJ' potential is supported.
        residues : str of list of str
            Labels of unique residues in the Compound. Residues are assigned by
            checking against Compound.name.  Only supply residue names as 3 character
            strings, as the residue names are truncated to 3 characters to fit in the
            psf and pdb file.
        forcefield_selection : str or dictionary, default = None
            Apply a forcefield to the output file by selecting a force field XML file with
            its path or by using the standard force field name provided the `foyer` package.
            Note: to write the NAMD/GOMC force field, pdb, and psf files, the
            residues and forcefields must be provided in a str or
            dictionary.  If a dictionary is provided all residues must
            be specified to a force field.
                Example dict for FF file: {'ETH' : 'oplsaa.xml', 'OCT': 'path_to_file/trappe-ua.xml'}
                Example str for FF file: 'path_to file/trappe-ua.xml'
                Example dict for standard FF names : {'ETH' : 'oplsaa', 'OCT': 'trappe-ua'}
                Example str for standard FF names: 'trappe-ua'
                Example of a mixed dict with both : {'ETH' : 'oplsaa', 'OCT': 'path_to_file/'trappe-ua.xml'}
        detect_forcefield_style: boolean
            If True, format NAMD/GOMC/LAMMPS parameters based on the contents of
            the parmed structure_box_0
        fix_res_bonds_angles: list, default = None
            When list of residues is provided, the selected residues will have
            their bonds and angles fixed in the GOMC engine.  This is specifically
            for the GOMC engine and it changes the residue's bond constants (Kbs)
            and angle constants (Kthetas) values to 999999999999 in the
            FF file (i.e., the .inp file).
        bead_to_atom_name_dict: dict, optional, default =None
            For all atom names/elements/beads with 2 or less digits, this converts
            the atom name in the GOMC psf and pdb files to a unique atom name,
            provided they do not exceed 3844 atoms (62^2) of the same name/element/bead
            per residue. For all atom names/elements/beads with 3 digits, this converts
            the atom name in the GOMC psf and pdb files to a unique atom name,
            provided they do not exceed 62 of the same name/element pre residue.
            Example dictionary: {'_CH3':'C', '_CH2':'C', '_CH':'C', '_HC':'C'}
        fix_residue: list  or None, default = None
            Changes occcur in the pdb file only.
            When residues are listed here, all the atoms in the residue are
            fixed and can not move via setting the Beta values in the PDB
            file to 1.00.
            If neither fix_residue or fix_residue_in_box lists a
            residue or both equal None, then the Beta values for all the atoms
            in the residue are free to move in the simulation and Beta values
            in the PDB file is set to 0.00
        fix_residue_in_box: list  or None, default = None
            Changes occcur in the pdb file only.
            When residues are listed here, all the atoms in the residue
            can move within the box but cannot be transferred between boxes
            via setting the Beta values in the PDB file to 2.00.
            If neither fix_residue or fix_residue_in_box lists a
            residue or both equal None, then the Beta values for all the atoms
            in the residue are free to move in the simulation and Beta values
            in the PDB file is set to 0.00
        FF_filename ; str, default =None
            If a string, it will write the  force field files that work in
            GOMC and NAMD
            structures
        reorder_res_in_pdb_psf ; bool, default =False
            If False, the order of of the atoms in the pdb file is kept in
            its original order, as in the Compound sent to the writer.
            If True, the order of the atoms is reordered based on their
            residue names in the 'residues' list that was entered.
        box_0 ; list of 3 positive float values or the dimensions [x, y ,z]
            for structure_box_0 in nanometers (nm)
            This is to add/override or change the structures dimensions. Ex: [1,2,3]
        box_1 ; list of 3 positive float values or the dimensions [x, y ,z]
            for structure_box_1 in nanometers (nm)
            This is to add/override or change the structures dimensions. Ex: [1,2,3]
        Notes
        -----
        Impropers and NBFIX are not currenly supported
        Currently the NBFIX is disabled as since only the OPLS and TRAPPE force fields are supported
        OPLS and CHARMM forcefield styles are supported, AMBER forcefield styles are NOT supported.
        Impropers and Urey-Bradleys are not supported for GOMC

        The atom typing is currently provided via a base 52 numbering (capital and lowercase lettering).
        This base 52 numbering allows for (52)^4 unique atom types.

        Unique atom names are provided if the system do not exceed 3844 atoms (62^2) of the same
        name/bead per residue (base 62 numbering). For all atom names/elements with 3 or less digits,
        this converts the atom name in the GOMC psf and pdb files to a unique atom name,
        provided they do not exceed 62 of the same name/element pre residue.

        Generating an empty pdb/psf:
            Single System: Enter residues = [], but the accompanying structure (structure_box_0)
            must be an empty mb.Compound structure. However, when doing this, the forcefield_selection
             must be supplied, or it will provide an error
            (i.e., forcefield_selection are both not equal to None)
            Dual System: Enter an empty mb.Compound structure for either structure_box_0 or
            structure_box_1.

        In this current FF/psf/pdb writer, a residue type is essentially a molecule type.
        Therefore, it can only correctly write systems where every bead/atom in the molecule
        has the same residue name, and the residue name is specific to that molecule type.
        For example: a protein molecule with many residue names is not currently supported,
        but is planned to be supported in the future.

        """

        # set all input variables to the class
        self.structure_box_0  = structure_box_0
        self.filename_box_0 = filename_box_0
        self.structure_box_1 = structure_box_1
        self.filename_box_1 = filename_box_1
        self.non_bonded_type =  non_bonded_type
        self.forcefield_selection = forcefield_selection
        self.residues = residues
        self.detect_forcefield_style = detect_forcefield_style
        self.fix_res_bonds_angles = fix_res_bonds_angles
        self.bead_to_atom_name_dict = bead_to_atom_name_dict
        self.fix_residue = fix_residue
        self.fix_residue_in_box = fix_residue_in_box
        self.FF_filename = FF_filename
        self.reorder_res_in_pdb_psf = reorder_res_in_pdb_psf
        self.box_0 = box_0
        self.box_1 = box_1
        #value to check for errors, with  self.input_error = True or False. Set to False initally
        self.input_error = False

        if not isinstance(self.structure_box_0, Compound):
            self.input_error = True
            print_error_message ='ERROR: The structure_box_0 expected to be of type: ' \
                                 '{}, received: {}'.format(type(Compound()), type(structure_box_0))
            raise TypeError(print_error_message)

        if  self.structure_box_1 != None and not isinstance(self.structure_box_1, Compound):
            self.input_error = True
            print_error_message = 'ERROR: The structure_box_1 expected to be of type: ' \
                                  '{}, received: {}'.format(type(Compound()), type(structure_box_1))
            raise TypeError(print_error_message)


        if not isinstance(self.residues, list):
            self.input_error = True
            print_error_message = 'ERROR: Please enter the residues list (residues) in a list format.'
            raise TypeError(print_error_message)

        if isinstance(self.residues, list):
            for each_residue in self.residues:
                if not isinstance(each_residue, str):
                    self.input_error = True
                    print_error_message = 'ERROR: Please enter a residues list (residues) with only string values.'
                    raise TypeError(print_error_message)

        if self.residues is None:
            self.input_error = True
            print_error_message = 'ERROR: Please enter the residues list (residues)'
            raise TypeError(print_error_message)
        if not isinstance(self.filename_box_0, str):
            self.input_error = True
            print_error_message = 'ERROR: Please enter the filename_box_0 as a string.'
            raise TypeError(print_error_message)

        unique_residue_test_name_list = []
        for res_m in range(0, len(self.residues)):
            if self.residues[res_m] not in unique_residue_test_name_list:
                unique_residue_test_name_list.append(self.residues[res_m])
        if len(unique_residue_test_name_list) != len(self.residues):
            self.input_error = True
            print_error_message = 'ERROR: Please enter the residues list (residues) that has only unique residue names.'
            raise ValueError(print_error_message)

        if self.filename_box_1 != None and not isinstance(self.filename_box_1, str):
            self.input_error = True
            print_error_message = 'ERROR: Please enter the filename_box_1 as a string.'
            raise TypeError(print_error_message)

        if self.structure_box_1 is None and self.box_1 != None:
            self.input_error = True
            print_error_message = 'ERROR: box_1 is set to a value but there is not a structure 1 to use it on.'
            raise TypeError(print_error_message)


        if self.FF_filename != None:
            if not isinstance(self.FF_filename, str):
                self.input_error = True
                print_error_message = 'ERROR: Please enter GOMC force field name (FF_filename) as a string.'
                raise TypeError(print_error_message)
            if isinstance(self.FF_filename, str):
                extension_FF_name = os.path.splitext(self.FF_filename)[-1]
                if extension_FF_name == '':
                    self.FF_filename = self.FF_filename + '.inp'
                elif extension_FF_name == '.inp':
                    self.FF_filename = self.FF_filename + ''
                elif extension_FF_name != '.inp':
                    self.input_error = True
                    print_error_message = 'ERROR: Please enter GOMC force field name without an '\
                                          'extention or the .inp extension.'
                    raise ValueError(print_error_message)

        if self.forcefield_selection != None:
            print('write_gomcdata: forcefield_selection = '+str(self.forcefield_selection) \
                  +', ' +'residues = '+str(self.residues) )
            if self.forcefield_selection != None and not isinstance(self.forcefield_selection, dict) \
                    and not isinstance(self.forcefield_selection, str):
                self.input_error = True
                print_error_message = 'ERROR: The force field selection (forcefield_selection) '\
                                      'is not a string or a dictionary with all the residues specified ' \
                                      'to a force field. -> String Ex: "path/trappe-ua.xml" or Ex: "trappe-ua" '\
                                      'Otherise provided a dictionary with all the residues specified ' \
                                      'to a force field '\
                                      '->Dictionary Ex: {"Water" : "oplsaa", "OCT": "path/trappe-ua.xml"}, ' \
                                      'Note: the file path must be specified the force field file if ' \
                                      'a standard foyer force field is not used.'
                raise TypeError(print_error_message)

            if isinstance(self.forcefield_selection, str):
                FF_name = self.forcefield_selection
                self.forcefield_selection = {}
                for i in range(0, len(self.residues)):
                    self.forcefield_selection.update({self.residues[i] : FF_name})
                print('FF forcefield_selection = '+str(self.forcefield_selection))

        elif self.forcefield_selection is None:
            self.input_error = True
            print_error_message = 'ERROR: Please enter the forcefield_selection as it was not provided.'
            raise TypeError(print_error_message)



        if self.residues != None and not isinstance(self.residues, list):
            self.input_error = True
            print_error_message = 'ERROR:  Please enter the residues (residues) in a list format'
            raise TypeError(print_error_message)


        if self.fix_res_bonds_angles != None and not isinstance(self.fix_res_bonds_angles, list):
            self.input_error = True
            print_error_message = 'ERROR: Please enter the residues that have fixed angles '\
                                  'and bonds (fix_res_bonds_angles) in a list format.'
            raise TypeError(print_error_message)

        if isinstance(self.fix_res_bonds_angles, list):
            for fix_res_bonds_angles_i in self.fix_res_bonds_angles:
                if fix_res_bonds_angles_i not in self.residues:
                    self.input_error = True
                    print_error_message = 'ERROR: Please ensure that all the residue names in the ' \
                                          'fix_res_bonds_angles list are also in the residues list.'
                    raise ValueError(print_error_message)
                elif not isinstance(fix_res_bonds_angles_i, str):
                    self.input_error = True
                    print_error_message = 'ERROR: Please enter a fix_res_bonds_angle list with only string values.'
                    raise TypeError(print_error_message)
                else:
                    print('INFORMATION: The following residues will have fixed bonds'
                          + ' and angles: fix_res_bonds_angles = ' +str(self.fix_res_bonds_angles))


        if self.fix_residue != None and not isinstance(self.fix_residue, list):
            self.input_error = True
            print_error_message = 'ERROR: Please enter the fix_residue in a list format'
            raise TypeError(print_error_message)

        if isinstance(self.fix_residue, list):
            for q in range(0,len(self.fix_residue)):
                if self.fix_residue[q] not in self.residues:
                    self.input_error = True
                    print_error_message = 'Error: Please ensure that all the residue names in the fix_residue '\
                                          'list are also in the residues list.'
                    raise ValueError(print_error_message)

        if self.fix_residue_in_box != None and not isinstance(self.fix_residue_in_box, list):
            self.input_error = True
            print_error_message = 'ERROR: Please enter the fix_residue_in_box in a list format.'
            raise TypeError(print_error_message)

        if isinstance(self.fix_residue_in_box, list):
            for q in range(0,len(self.fix_residue_in_box)):
                if self.fix_residue_in_box[q] not in self.residues:
                    self.input_error = True
                    print_error_message = 'Error: Please ensure that all the residue names in the ' \
                                          'fix_residue_in_box list are also in the residues list.'
                    raise ValueError(print_error_message)

        if self.bead_to_atom_name_dict != None and not isinstance(self.bead_to_atom_name_dict, dict):
            self.input_error = True
            print_error_message = 'ERROR: Please enter the a bead type to atom in the dictionary '\
                                  '(bead_to_atom_name_dict) so GOMC can properly evaluate the unique atom names'
            raise TypeError(print_error_message)

        if isinstance(self.bead_to_atom_name_dict, dict):
            dict_list = []
            for key in self.bead_to_atom_name_dict.keys():
                dict_list.append(key)

            for dict_lis_i in dict_list:
                if not isinstance(dict_lis_i, str) or not isinstance(self.bead_to_atom_name_dict[dict_lis_i], str):
                    print_error_message = 'ERROR: Please enter the bead_to_atom_name_dict with only string inputs.'
                    raise TypeError(print_error_message)

        if self.box_0 !=None :
            box_length = len(self.box_0)
            if box_length != 3:
                self.input_error = True
                print_error_message = 'ERROR: Please enter all 3 values and only 3 values for the box_0 dimensions.'
                raise ValueError(print_error_message)
            for box_iter in range(0, len(self.box_0)):
                for box_iter in range(0, len(self.box_0)):
                    if isinstance(self.box_0[box_iter], str) or self.box_0[box_iter] <= 0:
                        self.input_error = True
                        print_error_message = 'ERROR: Please enter float or integer values, which are all ' \
                                              'positive values for the box_0 dimensions.'
                        raise ValueError(print_error_message)

        if self.box_1 !=None :
            box_length = len(self.box_1)
            if box_length != 3:
                self.input_error = True
                print_error_message = 'ERROR: Please enter all 3 values and only 3 values for the box_1 dimensions.'
                raise ValueError(print_error_message)
            for box_iter in range(0, len(self.box_1)):
                if isinstance(self.box_1[ box_iter], str) or self.box_1[ box_iter] <= 0:
                    self.input_error = True
                    print_error_message = 'ERROR: Please enter float or integer values, which are all ' \
                                          'positive values for the box_1 dimensions.'
                    raise ValueError(print_error_message)


        print("******************************")
        print("")




        self.sub_1_for_base_52 = 1

        #if self.structure_box_1 != None:
        if self.structure_box_1:
            self.boxes_for_simulation = 2
        else:
            self.boxes_for_simulation = 1

        #write the Force fields
        self.combined_1_4_LJ_dict_per_residue = {}
        self.combined_1_4_Coul_dict_per_residue = {}
        #if self.structure_box_1 != None:
        if self.structure_box_1:

            print('GOMC FF writing each residues FF as a group for structure_box_0')
            self.structure_box_0_FF, \
            self.coulomb14scalar_dict_0, \
            self.LJ14scalar_dict_0,\
            self.residues_applied_list_0 = specific_ff_to_residue(self.structure_box_0,
                                                             forcefield_selection=self.forcefield_selection,
                                                             residues=self.residues,
                                                             reorder_res_in_pdb_psf=self.reorder_res_in_pdb_psf,
                                                             box = self.box_0,
                                                             boxes_for_simulation = self.boxes_for_simulation)
            test_Specific_FF_to_residue_for_failure = [self.structure_box_0_FF, self.coulomb14scalar_dict_0,
                                                       self.LJ14scalar_dict_0, self.residues_applied_list_0]

            for iter_test_Specifc_res_fail in range(0, len(test_Specific_FF_to_residue_for_failure)):
                if test_Specific_FF_to_residue_for_failure[iter_test_Specifc_res_fail] is None:
                    self.input_error = True
                    print_error_message = 'ERROR: The residues entered does not match the residues that were ' \
                                          'found and built for structure_box_0.'
                    raise ValueError(print_error_message)

            print('GOMC FF writing each residues FF as a group for  structure_box_1')
            self.structure_box_1_FF, \
            self.coulomb14scalar_dict_1, \
            self.LJ14scalar_dict_1, \
            self.residues_applied_list_1 = specific_ff_to_residue(self.structure_box_1,
                                                             forcefield_selection=self.forcefield_selection,
                                                             residues=self.residues,
                                                             reorder_res_in_pdb_psf=self.reorder_res_in_pdb_psf,
                                                             box = self.box_1,
                                                             boxes_for_simulation = self.boxes_for_simulation)
            test_Specific_FF_to_residue_for_failure = [self.structure_box_1_FF, self.coulomb14scalar_dict_1,
                                                       self.LJ14scalar_dict_1, self.residues_applied_list_1]

            for iter_test_Specifc_res_fail in range(0, len(test_Specific_FF_to_residue_for_failure)):
                if test_Specific_FF_to_residue_for_failure[iter_test_Specifc_res_fail] is None:
                    self.input_error = True
                    print_error_message = 'ERROR: The residues entered does not match the residues that were ' \
                                          'found and built for structure_box_1.'
                    raise ValueError(print_error_message)

            self.structure_box_0_and_1_FF =self.structure_box_0_FF + self.structure_box_1_FF
            self.combined_1_4_LJ_dict_per_residue.update(self.LJ14scalar_dict_0)
            self.combined_1_4_LJ_dict_per_residue.update(self.LJ14scalar_dict_1)
            self.combined_1_4_Coul_dict_per_residue.update(self.coulomb14scalar_dict_0)
            self.combined_1_4_Coul_dict_per_residue.update(self.coulomb14scalar_dict_1)

            self.residues_applied_list_0_and_1 = self.residues_applied_list_0
            for res_iter in range(0,len(self.residues_applied_list_1 )):
                if self.residues_applied_list_1[res_iter] not in self.residues_applied_list_0:
                    self.residues_applied_list_0_and_1.append(self.residues_applied_list_1[res_iter])

            for res_iter_0_1 in self.residues_applied_list_0_and_1:
                if res_iter_0_1 not in  self.residues:
                    self.input_error = True
                    print_error_message = "ERROR: All the residues were not used from the forcefield_selection " \
                                          "string or dictionary.  There may be residues below other specified " \
                                          "residues in the mbuild.Compound hierarchy.  If so, the residues " \
                                          "acquire the residue's force fields, which is at the top of the " \
                                          "hierarchy.  Alternatively, residues that are not in the structure " \
                                          "may have been specified."
                    raise ValueError(print_error_message)

            for res_iter_0_1 in self.residues:
                if res_iter_0_1 not in self.residues_applied_list_0_and_1:
                    self.input_error = True
                    print_error_message = ("ERROR: All the residues were not used from the forcefield_selection " \
                                           "string or dictionary.  There may be residues below other specified " \
                                           "residues in the mbuild.Compound hierarchy.  If so, the residues " \
                                           "acquire the residue's force fields, which is at the top of the " \
                                           "hierarchy.  Alternatively, residues that are not in the structure " \
                                           "may have been specified.")
                    raise ValueError(print_error_message)

            total_charge = sum([atom.charge for atom in self.structure_box_0_FF])
            if round(total_charge, 4) != 0.0:
                warn('System is not charge neutral for structure_box_0. '
                                 'Total charge is {}.'.format(total_charge))

            total_charge = sum([atom.charge for atom in self.structure_box_1_FF])
            if round(total_charge, 4) != 0.0:
                warn('System is not charge neutral for structure_box_1. '
                                 'Total charge is {}.'.format(total_charge))

            total_charge = sum([atom.charge for atom in self.structure_box_0_and_1_FF])
            if round(total_charge, 4) != 0.0:
                warn('System is not charge neutral for structure_0_and_1. '
                                 'Total charge is {}.'.format(total_charge))

        else:


            print('GOMC FF writing each residues FF as a group for structure_box_0')
            self.structure_box_0_FF, \
            self.coulomb14scalar_dict_0, \
            self.LJ14scalar_dict_0, \
            self.residues_applied_list_0 = specific_ff_to_residue(self.structure_box_0,
                                                             forcefield_selection=self.forcefield_selection,
                                                             residues=self.residues,
                                                             reorder_res_in_pdb_psf=self.reorder_res_in_pdb_psf,
                                                             box=self.box_0,
                                                             boxes_for_simulation = self.boxes_for_simulation)
            test_Specific_FF_to_residue_for_failure = [ self.structure_box_0_FF, self.coulomb14scalar_dict_0,
                                                        self.LJ14scalar_dict_0, self.residues_applied_list_0 ]
            for iter_test_Specifc_res_fail in range(0, len(test_Specific_FF_to_residue_for_failure)):
                if test_Specific_FF_to_residue_for_failure[iter_test_Specifc_res_fail] is None:
                    self.input_error = True
                    print_error_message = 'ERROR: The residues entered does not match the residues that were ' \
                                          'found and built for structure_box_0.'
                    raise ValueError(print_error_message)

            self.combined_1_4_LJ_dict_per_residue.update(self.LJ14scalar_dict_0)
            self.combined_1_4_Coul_dict_per_residue.update(self.coulomb14scalar_dict_0)

            for res_iter_0 in self.residues_applied_list_0:
                if res_iter_0 not in self.residues:
                    self.input_error = True
                    print_error_message = "ERROR: All the residues were not used from the forcefield_selection " \
                                           "string or dictionary.  There may be residues below other specified " \
                                           "residues in the mbuild.Compound hierarchy.  If so, the residues " \
                                           "acquire the residue's force fields, which is at the top of the " \
                                           "hierarchy.  Alternatively, residues that are not in the structure " \
                                           "may have been specified."
                    raise ValueError(print_error_message)

            for res_iter_0 in self.residues:
                if res_iter_0 not in self.residues_applied_list_0:
                    self.input_error = True
                    print_error_message = "ERROR: All the residues were not used from the forcefield_selection " \
                                           "string or dictionary.  There may be residues below other specified " \
                                           "residues in the mbuild.Compound hierarchy.  If so, the residues " \
                                           "acquire the residue's force fields, which is at the top of the " \
                                           "hierarchy.  Alternatively, residues that are not in the structure " \
                                           "may have been specified."
                    raise ValueError(print_error_message)

            total_charge = sum([atom.charge for atom in self.structure_box_0_FF])
            if round(total_charge, 4) != 0.0:
                warn('System is not charge neutral for structure_box_0. '
                                 'Total charge is {}.'.format(total_charge))


        print('forcefield type from compound = ' + str( self.forcefield_selection))
        print('coulomb14scale from compound = ' + str(self.combined_1_4_Coul_dict_per_residue))
        print('lj14scale from compound = ' + str(self.combined_1_4_LJ_dict_per_residue))

        # lock the atom_style and unit_style for GOMC. Can be inserted into variables once more functionality is built in
        self.atom_style = 'full'
        self.unit_style = 'real'
        # functional form type default.  Can be inserted into variables once more functionality is built in
        use_rb_torsions = True
        use_dihedrals = False
        use_urey_bradleys = False

        # Convert coordinates to real or other units (real only current option)
        if self.unit_style == 'real':
            self.sigma_conversion_factor = 1
            self.epsilon_conversion_factor = 1
            self.mass_conversion_factor = 1

        if self.structure_box_1:
            self.types = np.array([atom.type + '_' + str(atom.residue.name)
                                   for atom in self.structure_box_0_and_1_FF.atoms]
                                  )

        else:
            self.types = np.array([atom.type + '_' + str(atom.residue.name)
                                   for atom in self.structure_box_0_FF.atoms])

        self.unique_types = list(set(self.types))
        self.unique_types.sort(key=natural_sort)

        print('unique_types = {}'.format(self.unique_types))


        if self.structure_box_1:
            self.masses = np.array([atom.mass for atom in self.structure_box_0_and_1_FF.atoms]) / self.mass_conversion_factor
            self.mass_dict = dict([(self.unique_types.index(atom_type) + 1, mass) for atom_type, mass in zip(self.types, self.masses)])

        else:
            self.masses = np.array([atom.mass for atom in self.structure_box_0_FF.atoms]) / self.mass_conversion_factor
            self.mass_dict = dict([(self.unique_types.index(atom_type) + 1, mass) for atom_type, mass in zip(self.types, self.masses)])

        # added an index so the atom types can be converted to numbers as the type name is to long for insertion into
        # the pdb and psf files
        self.atom_types_to_index_value_dict = dict(
            [(self.unique_types[self.unique_types.index(atom_type)], self.unique_types.index(atom_type) )
             for atom_type, mass in zip(self.types, self.masses)])

        self.box_0 = Box(lengths=np.array([0.1 * val for val in self.structure_box_0_FF.box[0:3]]),angles=self.structure_box_0_FF.box[3:6])
        #Divide by conversion factor
        self.box_0.maxs /= self.sigma_conversion_factor

        #Internally use nm
        if self.structure_box_1:
            self.box_1 = Box(lengths=np.array([0.1 * val for val in self.structure_box_1_FF.box[0:3]]),
                      angles=self.structure_box_1_FF.box[3:6])
            # Divide by conversion factor
            self.box_1.maxs /= self.sigma_conversion_factor

        #if self.structure_box_1 != None:
        if self.structure_box_1:
            self.structure_selection = self.structure_box_0_and_1_FF
        else:
            self.structure_selection = self.structure_box_0_FF

        # Non-Bonded forces
        self.epsilons = np.array(
            [atom.epsilon for atom in self.structure_selection.atoms]) / self.epsilon_conversion_factor
        self.sigmas = np.array(
            [atom.sigma for atom in self.structure_selection.atoms]) / self.sigma_conversion_factor
        self.forcefields = [atom.type + '_' + atom.residue.name for atom in self.structure_selection.atoms]
        self.Residues = [atom.residue.name for atom in self.structure_selection.atoms]
        self.epsilon_dict = dict([(self.unique_types.index(atom_type), epsilon)
                                    for atom_type, epsilon in zip(self.types, self.epsilons)])
        self.sigma_dict = dict(
            [(self.unique_types.index(atom_type), sigma) for atom_type, sigma in zip(self.types, self.sigmas)])
        self.LJ_1_4_dict = dict(
            [(self.unique_types.index(atom_type), self.combined_1_4_LJ_dict_per_residue[self.Residues])
             for atom_type, self.Residues in zip(self.types, self.Residues)])
        self.forcefield_dict = dict([(self.unique_types.index(atom_type), forcefield)
                                       for atom_type, forcefield in zip(self.types, self.forcefields)])

        # ensure all 1,4-coulombic scaling factors are the same
        self.coul_1_4_List = []
        for p in self.combined_1_4_Coul_dict_per_residue.values():
            self.coul_1_4_List.append(p)
        self.coul_1_4_set = set(self.coul_1_4_List)
        if  len(self.coul_1_4_set) > 1:
            self.input_error = True
            print_error_message = "ERROR: There are multiple 1,4-coulombic scaling factors "\
                                  "GOMC will only accept a singular input for the 1,4-coulombic " \
                                  "scaling factors."
            raise ValueError(print_error_message)
        else:
            self.coul_1_4 = list(self.coul_1_4_set)[0]



        # get all the unique atom name to check for the MEMC move in the gomc_conf_writer
        self.all_Individual_atom_names_List = []
        self.all_residue_names_List = []
        if self.structure_box_1:
            list_of_structures = [self.structure_box_0_FF, self.structure_box_1_FF]
            stuct_only = [self.structure_box_0_FF, self.structure_box_1_FF]
        else:
            list_of_structures = [self.structure_box_0_FF]
            stuct_only = [self.structure_box_0_FF]

        for q_i in range(0, len(list_of_structures)):
            stuct_only_iteration =stuct_only[q_i]

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
                unique_residue_data_list.append(str(stuct_only_iteration.residues[m]))
                unique_residue_data_dict.update({unique_residue_data_list[m]: m + 1})
                residue_data_name_list.append(stuct_only_iteration.residues[m].name)

            self.Max_Residue_No = 9999
            self.No_1st_values_res_name = 3

            Res_No_iteration_corrected_List = []
            residue_ID_list = []
            for f, PSF_atom_iteration_0 in enumerate(stuct_only_iteration.atoms):

                residue_ID_int = int(unique_residue_data_dict[residue_data_list[f]])
                Res_ID_adder = int((residue_ID_int % self.Max_Residue_No) % (self.Max_Residue_No))
                if int(Res_ID_adder) == 0:
                    Res_No_iteration_corrected = int(self.Max_Residue_No)
                else:
                    Res_No_iteration_corrected = Res_ID_adder

                Res_No_iteration_corrected_List.append(Res_No_iteration_corrected)
                residue_ID_list.append(residue_ID_int)

            # This converts the atom name in the GOMC psf and pdb files to unique atom names
            unique_Individual_atom_names_dict_iter, \
            Individual_atom_names_List_iter, \
            Missing_Bead_to_atom_name_iter = unique_atom_naming(stuct_only_iteration ,
                                                               residue_ID_list,
                                                               residue_names_list,
                                                               bead_to_atom_name_dict=self.bead_to_atom_name_dict)

            print('Individual_atom_names_List_iter = {}'.format(Individual_atom_names_List_iter))
            print( 'self.all_Individual_atom_names_List = {}'.format(self.all_Individual_atom_names_List))
            if q_i == 0:
                self.all_Individual_atom_names_List = Individual_atom_names_List_iter

                self.all_residue_names_List = residue_names_list
            else:

                self.all_Individual_atom_names_List =  self.all_Individual_atom_names_List \
                                                              + Individual_atom_names_List_iter

                self.all_residue_names_List = self.all_residue_names_List \
                                               + residue_names_list

            print('Individual_atom_names_List_iter = {}'.format(Individual_atom_names_List_iter))
            print('self.all_Individual_atom_names_List = {}'.format(self.all_Individual_atom_names_List))

        # put the  self.all_Individual_atom_names_List and self.all_residue_names_List in a list to match
        # the the atom name with a residue and find the unique matches
        if None in [unique_Individual_atom_names_dict_iter, Individual_atom_names_List_iter, \
                    Missing_Bead_to_atom_name_iter]:
            self.input_error = True
            print_error_message = 'ERROR: The unique_atom_naming function failed while '\
                                  'running the charmm_writer function. Ensure the proper inputs are ' \
                                  'in the bead_to_atom_name_dict.'
            raise ValueError(print_error_message)

        else:
            self.all_atom_name_res_pairs_dict = {}
            for name_res_i in range(0, len(self.all_Individual_atom_names_List)):
                self.all_atom_name_res_pairs_dict.setdefault(self.all_residue_names_List[name_res_i], []
                                                             ).append(self.all_Individual_atom_names_List[name_res_i])

        print('self.all_atom_name_res_pairs_List = {}'.format(self.all_atom_name_res_pairs_dict))



    def write_inp(self):
        print("******************************")
        print("")
        print('The charmm force field file writer (the write_inp function) is running')

        if self.FF_filename is None:
            self.input_error = True
            print_error_message = 'ERROR: The force field file name was not specified and in the ' \
                                  'Charmm object. ' \
                                  'Therefore, the force field file (.inp) can not be written. ' \
                                  'Please use the force field file name when building the Charmm object, ' \
                                  'then use the write_inp function.'
            raise TypeError(print_error_message)
        else:

            print("******************************")
            print("")
            print('The charmm force field file writer (the write_inp function) is running')
            print("******************************")
            print("")
            print('writing the GOMC force field file ')
            date_time = datetime.datetime.today()

            unique_residue_data_dict = {}
            unique_residue_data_list = []
            residue_data_name_list = []
            #if self.structure_box_1 != None:
            if self.structure_box_1:
                residue_iterate = 0
                for m, residue in enumerate(self.structure_box_0_FF.residues):
                    residue_iterate = residue_iterate + 1
                    unique_residue_data_list.append(str(self.structure_box_0_FF.residues[m]))
                    unique_residue_data_dict.update({unique_residue_data_list[m]: m + 1})
                    residue_data_name_list.append(self.structure_box_0_FF.residues[m].name)

                for m, residue in enumerate(self.structure_box_1_FF.residues):
                    unique_residue_data_list.append(str(self.structure_box_1_FF.residues[m]))
                    unique_residue_data_dict.update({unique_residue_data_list[m]: m + 1 + residue_iterate})
                    residue_data_name_list.append(self.structure_box_1_FF.residues[m].name)


            else:
                for m, residue in enumerate(self.structure_box_0_FF.residues):
                    unique_residue_data_list.append(str(self.structure_box_0_FF.residues[m]))
                    unique_residue_data_dict.update({unique_residue_data_list[m]: m + 1})
                    residue_data_name_list.append(self.structure_box_0_FF.residues[m].name)

            for n in range(0, len(residue_data_name_list)):
                if residue_data_name_list[n] not in self.residues:
                    print('residue_data_name_list = ' + str(residue_data_name_list))
                    self.input_error = True
                    print_error_message = 'ERROR: Please specifiy all residues (residues) in a list'
                    raise ValueError(print_error_message)

            # Syntax which can change based on the functional form
            # Infer functional form based on the properties of the structure_box_0_and_1_FF or structure_box_0_FF
            if self.detect_forcefield_style:
                # Check angles
                if len(self.structure_selection.urey_bradleys) > 0 :
                    print("Urey bradley terms detected")
                    data.write("Urey bradley terms detected, will use angle_style self")
                    data.write("ERROR: GOMC does no support the Urey bradley terms")
                    use_urey_bradleys = True
                else:
                    print("No urey bradley terms detected, will use angle_style harmonic")
                    use_urey_bradleys = False

                # Check dihedrals
                if len(self.structure_selection.rb_torsions) > 0:
                    print("will use CHARMM_torsions  =  K0 + K1 * (1 + Cos[n1*(t) - (d1)] ) + "+
                          "K2 * (1 + Cos[n2*(t) - (d2)] ) + K3 * (1 + Cos[n3*(t) - (d3)] ) + "  +
                          "K4 * (1 + Cos[n4*(t) - (d4)] ) + K5 * (1 + Cos[n5*(t) - (d5)] ) ")
                    use_rb_torsions = True

                else:
                    use_rb_torsions = False

                if len(self.structure_selection.dihedrals) > 0:
                    print("Charmm dihedrals detected, will use dihedral_style charmm")
                    # this will need tested with a standard charmm input format before releasing it
                    use_dihedrals = True
                    self.input_error = True
                    print_error_message = "ERROR: use_dihedrals = {} " \
                                          "Charmm dihedrals not yet supported.".format(use_dihedrals)
                    raise ValueError(print_error_message)
                else:
                    use_dihedrals = False

            if use_rb_torsions and use_dihedrals:
                warn("Multiple dihedral styles detected, check your "
                                 "Forcefield XML and structure_selection")

            # Check impropers
            for dihedral in self.structure_selection.dihedrals:
                if dihedral.improper:
                    warn("Amber-style impropers are currently not supported")

            bonds = [[bond.atom1.idx+1, bond.atom2.idx+1] for bond in self.structure_selection.bonds]
            angles = [[angle.atom1.idx+1,
                       angle.atom2.idx+1,
                       angle.atom3.idx+1] for angle in self.structure_selection.angles]
            if use_rb_torsions:
                dihedrals = [[dihedral.atom1.idx+1,
                              dihedral.atom2.idx+1,
                              dihedral.atom3.idx+1,
                              dihedral.atom4.idx+1] for dihedral in self.structure_selection.rb_torsions]
            elif use_dihedrals:
                dihedrals = [[dihedral.atom1.idx+1,
                              dihedral.atom2.idx+1,
                              dihedral.atom3.idx+1,
                              dihedral.atom4.idx+1] for dihedral in self.structure_selection.dihedrals]
            else:
                dihedrals = []
            impropers = [[improper.atom1.idx+1,
                          improper.atom2.idx+1,
                          improper.atom3.idx+1,
                          improper.atom4.idx+1] for improper in self.structure_selection.impropers]


            if bonds :
                if len(self.structure_selection.bond_types) == 0:
                    bond_types = np.ones(len(bonds),dtype=int)
                else:
                    bond_types, unique_bond_types = _get_bond_types(self.structure_selection,
                                                                    self.sigma_conversion_factor,
                                                                    self.epsilon_conversion_factor)


            if angles:
                angle_types, unique_angle_types = _get_angle_types(
                    self.structure_selection,
                    self.sigma_conversion_factor,
                    self.epsilon_conversion_factor,
                    use_urey_bradleys=use_urey_bradleys
                )

            if dihedrals:
                dihedral_types, unique_dihedral_types = _get_dihedral_types(
                        self.structure_selection, use_rb_torsions, use_dihedrals,
                        self.epsilon_conversion_factor)



            if impropers:
                improper_types, unique_improper_types = _get_impropers(self.structure_selection,
                                                                       self.epsilon_conversion_factor)


            with open(self.FF_filename, 'w') as data:
                #if self.structure_box_1 != None:
                if self.structure_box_1:
                    data.write("*  " + self.filename_box_0 + ' and '+ self.filename_box_1+
                               ' - created by mBuild using the on ' + str(date_time) +'\n') #
                else:
                    data.write("*  " + self.filename_box_0 + ' - created by mBuild using the on '
                               + str(date_time) + '\n')

                data.write("*  " + 'parameters from the ' + str(self.forcefield_selection)
                           + ' force field(s) via MoSDef\n')
                data.write("*  1-4 coulombic scaling = " + str(self.combined_1_4_Coul_dict_per_residue)+
                           ', and 1-4 LJ scaling = ' + str(self.combined_1_4_LJ_dict_per_residue)+'\n\n')
                data.write("*  "+'{:d} atoms\n'.format(len(self.structure_selection.atoms)))

                if self.atom_style in ['full', 'molecular']:
                    data.write("*  "+'{:d} bonds\n'.format(len(bonds)))
                    data.write("*  "+'{:d} angles\n'.format(len(angles)))
                    data.write("*  "+'{:d} dihedrals\n'.format(len(dihedrals)))
                    data.write("*  "+'{:d} impropers\n\n'.format(len(impropers)))

                data.write("*  "+'{:d} atom types\n'.format(len(set(self.types))))
                if self.atom_style in ['full', 'molecular']:
                    if bonds:
                        data.write("*  "+'{:d} bond types\n'.format(len(set(unique_bond_types))))
                    if angles:
                        data.write("*  "+'{:d} angle types\n'.format(len(set(unique_angle_types))))
                    if dihedrals:
                        data.write("*  "+'{:d} dihedral types\n'.format(len(set(unique_dihedral_types))))
                    if impropers:
                        data.write("*  "+'{:d} improper types\n'.format(len(set(unique_improper_types))))


                data.write('\n')


                data.write('\n* masses\n\n')
                data.write('! atom_types \tmass \t\t  atomTypeForceFieldName_ResidueName '
                           +'(i.e., atoms_type_per_utilized_FF)  \n')
                for atom_type,mass in self.mass_dict.items():
                    mass_format = '*  {}\t\t{:.6f}\t! {}\n'
                    data.write(mass_format.format(base10_to_base52_alph(atom_type-self.sub_1_for_base_52),
                                                  mass, self.unique_types[atom_type - 1]))


                # Bond coefficients
                if bonds:
                    data.write('\n')
                    data.write('BONDS * harmonic\n')
                    data.write('!\n')
                    data.write('!V(bond) = Kb(b - b0)**2\n')
                    data.write('!\n')
                    data.write('!Kb: kcal/mole/A**2\n')
                    data.write('!b0: A\n')
                    data.write('!Kb (kcal/mol) = Kb_K (K) * Boltz. const.; (9999999999 if no stretching)\n')
                    data.write('!\n')


                    if self.unit_style == 'real':
                        data.write('!atom_types \t Kb\t\tb0 \t\t  atoms_types_per_utilized_FF\n')
                    for params, idx in unique_bond_types.items():
                        bond_format = '{}\t{}\t{}\t{}\t\t! {}\t{}\n'
                        if (self.fix_res_bonds_angles != None) and ((params[3] and  params[4])
                                                                      in self.fix_res_bonds_angles ):
                            fix_bond_K_value = '999999999999'
                            data.write( bond_format.format(base10_to_base52_alph(self.atom_types_to_index_value_dict[params[2][0]+'_' + str(params[3])]),
                                                           base10_to_base52_alph(self.atom_types_to_index_value_dict[params[2][1]+'_' + str(params[4])]),
                                                           fix_bond_K_value, params[1],
                                                           params[2][0]+'_' + str(params[3]),
                                                           params[2][1]+'_' + str(params[4])))

                        else:
                            data.write( bond_format.format(base10_to_base52_alph(self.atom_types_to_index_value_dict[params[2][0] +'_' + str(params[3])]),
                                                           base10_to_base52_alph(self.atom_types_to_index_value_dict[params[2][1]+'_' + str(params[4])]),
                                                           params[0],params[1],
                                                           params[2][0]+'_' + str(params[3]),
                                                           params[2][1]+'_' + str(params[4])))


                # Angle coefficients
                if angles:
                    if use_urey_bradleys:
                        data.write('\n!  Urey Bradley terms detected but not written,'+
                                   'since they are currently not compatible with GOMC\n')

                    data.write('\nANGLES * harmonic\n')
                    data.write('!\n')
                    data.write('!V(angle) = Ktheta(Theta - Theta0)**2\n')
                    data.write('!\n')
                    data.write('!Ktheta: kcal/mole/rad**2\n')
                    data.write('!Theta0: degrees\n')
                    data.write('!\n')
                    data.write('! Ktheta (kcal/mol) = Ktheta_K (K) * Boltz. const.\t\t\n')
                    data.write('!\n')
                    data.write('!atom_types \t\tKtheta\t\tTheta0\t\t\t  atoms_types_per_utilized_FF\n')
                    for params,idx in unique_angle_types.items():

                        if (self.fix_res_bonds_angles != None) and ((params[4] and  params[5] and  params[6])
                                                               in self.fix_res_bonds_angles ):
                            fix_angle_K_value = '999999999999'
                            angle_format = '{}\t{}\t{}\t{}\t\t{:.5f}\t\t! {}\t{}\t{}\n'
                            data.write(angle_format.format(base10_to_base52_alph(self.atom_types_to_index_value_dict[params[3][0]+'_'+params[4]]),
                                                           base10_to_base52_alph(self.atom_types_to_index_value_dict[params[2]+'_'+params[5]]),
                                                           base10_to_base52_alph(self.atom_types_to_index_value_dict[params[3][1]+'_'+params[6]]),
                                                           fix_angle_K_value ,params[1],
                                                           params[3][0]+'_'+params[4],
                                                           params[2]+'_'+params[5],
                                                           params[3][1]+'_'+params[6]))

                        else:
                            data.write(
                                '{}\t{}\t{}\t{}\t\t{:.5f}\t\t! {}\t{}\t{}\n'.format(base10_to_base52_alph(self.atom_types_to_index_value_dict[params[3][0]+'_'+params[4]]),
                                                                                  base10_to_base52_alph(self.atom_types_to_index_value_dict[params[2]+'_'+params[5]]),
                                                                                  base10_to_base52_alph(self.atom_types_to_index_value_dict[params[3][1]+'_'+params[6]]),
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
                        data.write('! Kchi (kcal/mol) = Kchi_K (K) * Boltz. const.\n')
                        data.write('! Boltzmann = 0.0019872041 kcal / (mol * K)\n')
                        data.write('!\n')
                        if self.unit_style == 'real':
                            data.write('!atom_types \t\t\tKchi\t\tn\tdelta\t\t  atoms_types_per_utilized_FF\n')

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
                            data.write(dihedral_format.format(base10_to_base52_alph(self.atom_types_to_index_value_dict[ params[8]+ '_' + params[12]]),
                                                              base10_to_base52_alph(self.atom_types_to_index_value_dict[params[9]+ '_' + params[13]]),
                                                              base10_to_base52_alph(self.atom_types_to_index_value_dict[params[10]+ '_' + params[14]]),
                                                              base10_to_base52_alph(self.atom_types_to_index_value_dict[params[11]+ '_' + params[15]]),
                                                              CHARMM_coeffs[0,0],
                                                              int(CHARMM_coeffs[0,1]),
                                                              CHARMM_coeffs[0,2],
                                                              params[8]+ '_' + params[12],
                                                              params[9]+ '_' + params[13],
                                                              params[10]+ '_' + params[14],
                                                              params[11]+ '_' + params[15]))
                            data.write(dihedral_format.format(base10_to_base52_alph(self.atom_types_to_index_value_dict[params[8]+ '_' + params[12]]),
                                                              base10_to_base52_alph(self.atom_types_to_index_value_dict[params[9]+ '_' + params[13]]),
                                                              base10_to_base52_alph(self.atom_types_to_index_value_dict[params[10]+ '_' + params[14]]),
                                                              base10_to_base52_alph(self.atom_types_to_index_value_dict[params[11]+ '_' + params[15]]),
                                                              CHARMM_coeffs[1, 0],
                                                              int(CHARMM_coeffs[1, 1]),
                                                              CHARMM_coeffs[1, 2],
                                                              params[8]+ '_' + params[12],
                                                              params[9]+ '_' + params[13],
                                                              params[10]+ '_' + params[14],
                                                              params[11]+ '_' + params[15]))
                            data.write(dihedral_format.format(base10_to_base52_alph(self.atom_types_to_index_value_dict[params[8]+ '_' + params[12]]),
                                                              base10_to_base52_alph(self.atom_types_to_index_value_dict[params[9]+ '_' + params[13]]),
                                                              base10_to_base52_alph(self.atom_types_to_index_value_dict[params[10]+ '_' + params[14]]),
                                                              base10_to_base52_alph(self.atom_types_to_index_value_dict[params[11]+ '_' + params[15]]),
                                                              CHARMM_coeffs[2, 0],
                                                              int(CHARMM_coeffs[2, 1]),
                                                              CHARMM_coeffs[2, 2],
                                                              params[8]+ '_' + params[12],
                                                              params[9]+ '_' + params[13],
                                                              params[10]+ '_' + params[14],
                                                              params[11]+ '_' + params[15]))
                            data.write(dihedral_format.format(base10_to_base52_alph(self.atom_types_to_index_value_dict[params[8]+ '_' + params[12]]),
                                                              base10_to_base52_alph(self.atom_types_to_index_value_dict[params[9]+ '_' + params[13]]),
                                                              base10_to_base52_alph(self.atom_types_to_index_value_dict[params[10]+ '_' + params[14]]),
                                                              base10_to_base52_alph(self.atom_types_to_index_value_dict[params[11]+ '_' + params[15]]),
                                                              CHARMM_coeffs[3, 0],
                                                              int(CHARMM_coeffs[3, 1]),
                                                              CHARMM_coeffs[3, 2],
                                                              params[8]+ '_' + params[12],
                                                              params[9]+ '_' + params[13],
                                                              params[10]+ '_' + params[14],
                                                              params[11]+ '_' + params[15]))
                            data.write(dihedral_format.format(base10_to_base52_alph(self.atom_types_to_index_value_dict[params[8]+ '_' + params[12]]),
                                                              base10_to_base52_alph(self.atom_types_to_index_value_dict[params[9]+ '_' + params[13]]),
                                                              base10_to_base52_alph(self.atom_types_to_index_value_dict[params[10]+ '_' + params[14]]),
                                                              base10_to_base52_alph(self.atom_types_to_index_value_dict[params[11]+ '_' + params[15]]),
                                                              CHARMM_coeffs[4, 0],
                                                              int(CHARMM_coeffs[4, 1]),
                                                              CHARMM_coeffs[4, 2],
                                                              params[8]+ '_' + params[12],
                                                              params[9]+ '_' + params[13],
                                                              params[10]+ '_' + params[14],
                                                              params[11]+ '_' + params[15]))
                            data.write(dihedral_format.format(base10_to_base52_alph(self.atom_types_to_index_value_dict[params[8]+ '_' + params[12]]),
                                                              base10_to_base52_alph(self.atom_types_to_index_value_dict[params[9]+ '_' + params[13]]),
                                                              base10_to_base52_alph(self.atom_types_to_index_value_dict[params[10]+ '_' + params[14]]),
                                                              base10_to_base52_alph(self.atom_types_to_index_value_dict[params[11]+ '_' + params[15]]),
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
                            info_if_dihedral_error_too_large = '! WARNING: RB-torsion to CHARMM ' \
                                                               'dihedral conversion error' \
                                                               ' is to large [error > 10^(-10)] \n' \
                                                               '! WARNING: Maximum( ' \
                                                               '|(RB-torsion calc)-(CHARMM dihedral calc)| ) =  ' \
                                                               + str(List_if_largest_error_abs_Values_for_dihedral_overall_max)  \
                                                               +'\n'
                            warn('! WARNING: RB-torsion to CHARMM dihedral conversion error'
                                 'is to large [error > 10^(-10)] \n'
                                 '! WARNING: Maximum( |(RB-torsion calc)-(CHARMM dihedral calc)| ) =  ' \
                                 + str( List_if_largest_error_abs_Values_for_dihedral_overall_max)  + '\n')
                            data.write(info_if_dihedral_error_too_large)
                            print(info_if_dihedral_error_too_large)
                        else:
                            List_if_abs_max_Values_for_dihedral_overall_max = \
                                max(List_if_abs_max_Values_for_dihedral_overall)
                            info_if_dihedral_error_OK = '! RB-torsion to CHARMM dihedral conversion error is OK '\
                                                        '[error <= 10^(-10)]\n' + \
                                                        '! Maximum( |(RB-torsion calc)-(CHARMM dihedral calc)| ) =  ' \
                                                        + str(List_if_abs_max_Values_for_dihedral_overall_max) + '\n'
                            data.write(info_if_dihedral_error_OK)
                            print(info_if_dihedral_error_OK)




                    elif use_dihedrals:
                        data.write("ERROR: not set up to use to use_dihedrals form for data input from the xml file")

                # Improper coefficients
                if impropers:
                    data.write("ERROR: GOMC is not currently able to use improper in its calculations")


                # Pair coefficients
                print('NBFIX_Mixing not used or no mixing used for the non-bonded potentials out')
                print('self.non_bonded_type = ' +str(self.non_bonded_type))
                if self.non_bonded_type=='LJ':
                    data.write('\n')
                    data.write('NONBONDED\n')
                    data.write('!\n')
                    data.write('!V(Lennard-Jones) = Eps,i,j[(Rmin,i,j/ri,j)**12 - 2(Rmin,i,j/ri,j)**6]\n')
                    data.write('!\n')

                    if self.unit_style == 'real':
                        data.write('!atype \tignored\tepsilon \tRmin/2 \t\tignored\teps,1-4\t\tRmin/2,1-4\t\t'+
                                   '  atom_type_per_utilized_FF\n')

                    print('forcefield_dict = '+str(self.forcefield_dict))

                    for idx, epsilon in self.epsilon_dict.items():
                        NB_format = '{}\t{:.2f}\t{:.9f}\t{:.11f}\t{:.2f}\t{:.9f}\t{:.11f}\t\t! {}\t{}\n'
                        data.write(NB_format.format(base10_to_base52_alph(idx), 0, -epsilon,
                                                    self.sigma_dict[idx] * (2 ** (1 / 6)) / 2, 0,
                                                    float(self.LJ_1_4_dict[idx])* (-epsilon),
                                                    self.sigma_dict[idx] * (2 ** (1 / 6)) / 2,
                                                    self.forcefield_dict[idx],self.forcefield_dict[idx]))


                elif  self.non_bonded_type=='Mie':
                    data.write("ERROR: Currenly the Mie potential (non_bonded_type) is not supported in this MoSDeF "
                               "GOMC parameter writer\n")
                    print_error_message = "ERROR: Currenly the Mie potential (non_bonded_type) is not " \
                                          "supported in this MoSDeF GOMC parameter writer."
                    raise ValueError(print_error_message)
                else:
                    data.write("ERROR: Currenly this potential (non_bonded_type) is not supported in "
                               "this MoSDeF GOMC parameter writer\n")
                    print_error_message = "ERROR: Currenly this potential (non_bonded_type) is not supported in " \
                                          "this MoSDeF GOMC parameter writer."
                    raise ValueError(print_error_message)

                # writing end in file
                data.write('\nEND\n')



        # **********************************
        #**********************************
        # FF writer (end)
        # **********************************
        # **********************************


    def write_psf(self):
        # **********************************
        #**********************************
        # psf writer (start)
        # **********************************
        # **********************************


        print("******************************")
        print("")
        print('The charmm X-plor format psf writer (the write_psf function) is running')

        date_time = datetime.datetime.today()


        print('write_psf: forcefield_selection = {}, residues = {}'.format(self.forcefield_selection, self.residues))



        print("******************************")
        print("")


        if self.structure_box_1:
            list_of_structures = [self.structure_box_0_FF, self.structure_box_1_FF]
            list_of_file_names = [self.filename_box_0, self.filename_box_1]
            stuct_only = [self.structure_box_0_FF, self.structure_box_1_FF]
        else:
            list_of_structures = [self.structure_box_0_FF]
            list_of_file_names = [self.filename_box_0]
            stuct_only = [self.structure_box_0_FF]

        for q in range(0, len(list_of_structures)):
            stuct_iteration = list_of_structures[q]
            file_name_iteration = list_of_file_names[q]
            output = str(file_name_iteration)+'.psf'
            stuct_only_iteration =stuct_only[q]
            # Lammps syntax depends on the functional form
            # Infer functional form based on the properties of the stuct_iteration
            if self.detect_forcefield_style:
                # Check  for angles
                if len(stuct_iteration.urey_bradleys) > 0:
                    print("Warning: Urey bradley terms detected. GOMC does no support the Urey-Bradley terms")
                    warn("warning: Urey bradley terms detected. "
                                     "GOMC does no support the Urey-Bradley terms")
                    use_urey_bradleys = True
                else:
                    print("No urey bradley terms detected")
                    use_urey_bradleys = False

                # Check for dihedrals
                if len(stuct_iteration.rb_torsions) > 0:
                    print("RB Torsions detected, will converted to CHARMM Dihedrals")
                    use_rb_torsions = True
                    Dihedrals_list = stuct_iteration.rb_torsions
                    dihedrals = [[dihedral.atom1.idx + 1,
                                  dihedral.atom2.idx + 1,
                                  dihedral.atom3.idx + 1,
                                  dihedral.atom4.idx + 1] for dihedral in stuct_iteration.rb_torsions]
                else:
                    use_rb_torsions = False

                if len(stuct_iteration.dihedrals) > 0:
                    print("Charmm dihedrals detected, so CHARMM Dihedrals will remain")
                    use_dihedrals = True
                    Dihedrals_list = stuct_iteration.dihedrals
                    dihedrals = [[dihedral.atom1.idx + 1,
                                  dihedral.atom2.idx + 1,
                                  dihedral.atom3.idx + 1,
                                  dihedral.atom4.idx + 1] for dihedral in stuct_iteration.dihedrals]
                else:
                    use_dihedrals = False
            if (use_rb_torsions == False) and (use_dihedrals == False):
                Dihedrals_list = []
                dihedrals = []
            if use_rb_torsions and use_dihedrals:
                warn("Multiple dihedral styles detected, check your "
                                 "Forcefield XML and structure files")

            # Check for impropers
            for dihedral in stuct_iteration.dihedrals:
                if dihedral.improper:
                    warn("ERROR: Amber-style impropers are currently not supported in GOMC")

            impropers = [[improper.atom1.idx + 1,
                          improper.atom2.idx + 1,
                          improper.atom3.idx + 1,
                          improper.atom4.idx + 1] for improper in stuct_iteration.impropers]

            No_atoms = len(stuct_iteration.atoms)
            No_bonds = len(stuct_iteration.bonds)
            No_angles = len(stuct_iteration.angles)

            No_dihedrals = len(dihedrals)
            No_impropers = len(impropers)

            No_donors = len(stuct_iteration.donors)
            No_acceptors = len(stuct_iteration.acceptors)
            No_groups = len(stuct_iteration.groups)


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
                unique_residue_data_list.append(str(stuct_only_iteration.residues[m]))
                unique_residue_data_dict.update({unique_residue_data_list[m]: m + 1})
                residue_data_name_list.append(stuct_only_iteration.residues[m].name)


            Res_No_iteration_corrected_List = []
            residue_ID_list = []
            for f, PSF_atom_iteration_0 in enumerate(stuct_only_iteration.atoms):

                residue_ID_int = int(unique_residue_data_dict[residue_data_list[f]])
                Res_ID_adder = int((residue_ID_int % self.Max_Residue_No) % (self.Max_Residue_No))
                if int(Res_ID_adder) == 0:
                    Res_No_iteration_corrected = int(self.Max_Residue_No)
                else:
                    Res_No_iteration_corrected = Res_ID_adder

                Res_No_iteration_corrected_List.append(Res_No_iteration_corrected)
                residue_ID_list.append(residue_ID_int)


            output_write = genopen(output, 'w')

            first_indent = '%8s'
            PSF_formating = ('%8s %-4s %-4s %-4s %-4s %4s %10.6f %13.4f' + 11 * ' ')

            output_write.write('PSF ')
            output_write.write('\n\n')

            No_of_remarks = 3
            output_write.write(first_indent % No_of_remarks + ' !NTITLE\n')
            output_write.write(' REMARKS this file ' + file_name_iteration
                               + ' - created by mBuild/foyer using the' + '\n')
            output_write.write(' REMARKS parameters from the ' + str(self.forcefield_selection)
                               + ' force field via MoSDef\n')
            output_write.write(' REMARKS created on ' + str(date_time) + '\n\n\n')

            # This converts the atom name in the GOMC psf and pdb files to unique atom names
            print('bead_to_atom_name_dict = {}'.format(self.bead_to_atom_name_dict))
            unique_Individual_atom_names_dict, \
            Individual_atom_names_List, \
            Missing_Bead_to_atom_name = unique_atom_naming(stuct_only_iteration , residue_ID_list, residue_names_list,
                                                          bead_to_atom_name_dict=self.bead_to_atom_name_dict)

            if None in [unique_Individual_atom_names_dict, Individual_atom_names_List, Missing_Bead_to_atom_name]:
                self.input_error = True
                print_error_message = 'ERROR: The unique_atom_naming function failed while '\
                                      'running the charmm_writer function. Ensure the proper inputs are ' \
                                      'in the bead_to_atom_name_dict.'
                raise ValueError(print_error_message)

            # ATOMS: Calculate the atom data
            # PSF_formating is conducted for the for CHARMM format (i.e., atom types are base 52, letters only)
            output_write.write(first_indent % No_atoms  + ' !NATOM\n')
            for i_atom, PSF_atom_iteration_1 in enumerate(stuct_iteration.atoms):
                Segment_ID = PSF_atom_iteration_1.residue.segid or 'SYS'
                atom_type_iter = base10_to_base52_alph(self.atom_types_to_index_value_dict[PSF_atom_iteration_1.type
                                                                                           + '_' +
                                                                                          PSF_atom_iteration_1.residue.name])

                atom_lines_iteration = PSF_formating % (i_atom + 1, Segment_ID, Res_No_iteration_corrected_List[i_atom],
                                                        str(residue_names_list[i_atom])[:self.No_1st_values_res_name],
                                                        Individual_atom_names_List[i_atom], atom_type_iter,
                                                        PSF_atom_iteration_1.charge, PSF_atom_iteration_1.mass)

                output_write.write('%s\n' % atom_lines_iteration)

            output_write.write('\n')

            # BONDS: Calculate the bonding data
            output_write.write(first_indent % No_bonds + ' !NBOND: bonds\n')
            for i_bond, PSF_bond_iteration_1 in enumerate(stuct_iteration.bonds):
                output_write.write((first_indent * 2) % (PSF_bond_iteration_1.atom1.idx + 1,
                                                         PSF_bond_iteration_1.atom2.idx + 1))

                if (i_bond + 1) % 4 == 0:
                    output_write.write('\n')

            if No_bonds % 4 == 0:
                output_write.write('\n')
            else:
                output_write.write('\n\n')

            if No_bonds == 0:
                output_write.write('\n')

            # ANGLES: Calculate the angle data
            output_write.write(first_indent % No_angles + ' !NTHETA: angles\n')
            for i_angle, angle_iteration in enumerate(stuct_iteration.angles):

                output_write.write((first_indent * 3) % (angle_iteration.atom1.idx + 1, angle_iteration.atom2.idx + 1,
                                                         angle_iteration.atom3.idx + 1) )

                if (i_angle + 1) % 3 == 0:
                    output_write.write('\n')

            if No_angles % 3 == 0:
                output_write.write('\n')
            else:
                output_write.write('\n\n')

            if No_angles == 0:
                output_write.write('\n')

            # DIHEDRALS: Calculate the dihedral  data
            output_write.write(first_indent % No_dihedrals + ' !NPHI: dihedrals\n')
            for i_dihedral, dihedral_iter in enumerate(Dihedrals_list):
                dihedral_atom_1,  dihedral_atom_2,  dihedral_atom_3,  dihedral_atom_4 =  dihedral_iter.atom1,  \
                                                                                         dihedral_iter.atom2,  \
                                                                                         dihedral_iter.atom3,  \
                                                                                         dihedral_iter.atom4

                output_write.write((first_indent * 4) % (dihedral_atom_1.idx + 1, dihedral_atom_2.idx + 1,
                                                         dihedral_atom_3.idx+ 1, dihedral_atom_4.idx + 1))

                if (i_dihedral +1) % 2 == 0:
                    output_write.write('\n')

            if No_dihedrals % 2 == 0:
                output_write.write('\n')
            else:
                output_write.write('\n\n')

            if No_dihedrals == 0:
                output_write.write('\n')

            # IMPROPERS: Calculate the improper data
            output_write.write(first_indent % No_impropers + ' !NIMPHI: impropers\n')
            for i_improper, (atom_1, atom_2, atom_3, atom_4) in enumerate(generate_impropers_for_PSF(stuct_iteration,
                                                                                                     Dihedrals_list )):

                output_write.write((first_indent * 4) % (atom_1.idx + 1, atom_2.idx + 1, atom_3.idx + 1, atom_4.idx + 1))
                if (i_improper +1) % 2 == 0:
                    output_write.write('\n')

            if No_impropers % 2 == 0:
                output_write.write('\n')
            else:
                output_write.write('\n\n')

            if No_impropers == 0:
                output_write.write('\n')

            # DONOR: calculate the donor data
            output_write.write(first_indent % No_donors + ' !NDON: donors\n')
            for donor_i, donor_iter in enumerate(stuct_iteration.donors):

                output_write.write((first_indent * 2) % (donor_iter.atom1.idx + 1, donor_iter.atom2.idx + 1))
                if (donor_i + 1) % 4 == 0:
                    output_write.write('\n')

            if No_donors % 4 == 0:
                output_write.write('\n')
            else:
                output_write.write('\n\n')

            if No_donors == 0:
                output_write.write('\n')

            # ACCEPTOR: calculate the acceptor data
            output_write.write(first_indent % No_acceptors + ' !NACC: acceptors\n')
            for acceptor_i, acceptor_iter in enumerate(stuct_iteration.acceptors):

                output_write.write((first_indent * 2) % (acceptor_iter.atom1.idx + 1, acceptor_iter.atom2.idx + 1))
                if (acceptor_i + 1) % 4 == 0:
                    output_write.write('\n')

            if No_acceptors % 4 == 0:
                output_write.write('\n')
            else:
                output_write.write('\n\n')

            if No_acceptors == 0:
                output_write.write('\n')

            # NNB: calculate the NNB data
            output_write.write(first_indent % 0 + ' !NNB\n\n')
            for nbb_i, atoms_iter in enumerate(stuct_iteration.atoms):

                output_write.write(first_indent % 0)
                if (nbb_i + 1) % 8 == 0:
                    output_write.write('\n')

            if No_atoms % 8 == 0:
                output_write.write('\n')
            else:
                output_write.write('\n\n')

            if No_atoms == 0:
                output_write.write('\n')

            # GROUP: calculate the group data
            try:
                group_data = stuct_iteration.groups.nst2
            except AttributeError:
                group_data = 0
            output_write.write((first_indent * 2) % (No_groups or 1, group_data) + ' !NGRP \n')
            if stuct_iteration.groups == True:
                for group_i, group_iter in enumerate(stuct_iteration.groups):

                    output_write.write((first_indent * 3) % (group_iter.atom.idx, group_iter.type, group_iter.move))
                    if (group_i + 1) % 3 == 0:
                        output_write.write('\n')

                if No_groups % 3 == 0:
                    output_write.write('\n')
                else:
                    output_write.write('\n\n')

                if No_groups == 0:
                    output_write.write('\n')

            else:
                structure_abs_charge_value = abs(sum(atom_charge_iter.charge for atom_charge_iter in stuct_iteration.atoms))
                if structure_abs_charge_value < 1.0e-4:
                    group_type = 1
                else:
                    group_type = 2
                output_write.write((first_indent * 3) % (0, group_type, 0))
                output_write.write('\n')

            output_write.write('\n')
            output_write.close()
        # **********************************
        # **********************************
        # psf writer (end)
        # **********************************
        # **********************************


    def write_pdb(self):
        # **********************************
        # **********************************
        # pdb writer (start)
        # **********************************
        # **********************************
        date_time = datetime.datetime.today()
        print("******************************")
        print("")
        print('The charmm pdb writer (the write_pdb function) is running')
        print('write_charmm_pdb: residues == {}'.format(self.residues))
        print('fix_residue = {}'.format(self.fix_residue))
        print('fix_residue_in_box = {}'.format(self.fix_residue_in_box))
        print('bead_to_atom_name_dict = {}'.format(self.bead_to_atom_name_dict))

        if self.fix_residue is None and self.fix_residue_in_box is None:
            print('INFORMATION: No atoms are fixed in this pdb file for the GOMC simulation engine. ')
        else:
            warn('Some atoms are fixed in this pdb file for the GOMC simulation engine. ')

        print("******************************")
        print("")

        if self.structure_box_1:
            list_of_structures = [self.structure_box_0_FF, self.structure_box_1_FF]
            list_of_file_names = [self.filename_box_0, self.filename_box_1]
            stuct_only = [self.structure_box_0_FF, self.structure_box_1_FF]
        else:
            list_of_structures = [self.structure_box_0_FF]
            list_of_file_names = [self.filename_box_0]
            stuct_only = [self.structure_box_0_FF]

        for q in range(0, len(list_of_structures)):
            file_name_iteration = list_of_file_names[q]
            output = str(file_name_iteration)+'.pdb'
            stuct_only_iteration =stuct_only[q]

            output_write = genopen(output, 'w')
            output_write.write(
                'REMARK this file ' + file_name_iteration + ' - created by mBuild/foyer using the' + '\n')
            output_write.write(
                'REMARK parameters from the ' + str(self.forcefield_selection) + ' force field via MoSDef\n')
            output_write.write('REMARK created on ' + str(date_time) + '\n')

            unique_residue_data_dict = {}
            unique_residue_data_list = []
            residue_data_name_list = []
            for m, residue in enumerate( stuct_only_iteration.residues):
                unique_residue_data_list.append(str( stuct_only_iteration.residues[m]))
                unique_residue_data_dict.update({ unique_residue_data_list[m]: m+1})
                residue_data_name_list.append( stuct_only_iteration.residues[m].name)

            for n in range(0, len(residue_data_name_list)):
                if residue_data_name_list[n] not in  self.residues:
                    self.input_error = True
                    print_error_message = 'ERROR: Please specifiy all residues (residues) in a list'
                    raise ValueError(print_error_message)

            residue_data_list = []
            for k, atom in enumerate( stuct_only_iteration.atoms):
                residue_data_list.append(str(atom.residue))

            if (self.fix_residue != None) and (self.fix_residue_in_box != None):
                for n in range(0,len(self.fix_residue)):
                    if self.fix_residue[n] in self.fix_residue_in_box:
                        self.input_error = True
                        print_error_message = "ERROR: residue type can not be specified to both "\
                                              "fix_residue and fix_residue_in_box"
                        raise ValueError(print_error_message)

            residue_names_list =[]
            fix_atoms_list = []
            for k, atom in enumerate( stuct_only_iteration.atoms):
                residue_names_list.append(atom.residue.name)
                if (self.fix_residue != None) and (atom.residue.name  in self.fix_residue) :
                        beta_iteration = 1.00
                elif (self.fix_residue_in_box != None) and (atom.residue.name  in self.fix_residue_in_box):
                    beta_iteration = 2.00
                else:
                    beta_iteration = 0.00
                fix_atoms_list.append(beta_iteration)

            if  stuct_only_iteration.box is not None:
                output_write.write('CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f %-11s%4s\n' % (stuct_only_iteration.box[0],
                                                                                        stuct_only_iteration.box[1],
                                                                                        stuct_only_iteration.box[2],
                                                                                        stuct_only_iteration.box[3],
                                                                                        stuct_only_iteration.box[4],
                                                                                        stuct_only_iteration.box[5],
                                                                                        stuct_only_iteration.space_group,
                                                                                  ''
                                                                                  ))

            All_atom_coordinates = stuct_only_iteration.get_coordinates('all')
            if All_atom_coordinates is None:
                self.input_error = True
                print_error_message = "ERROR: the submitted structure has no PDB coordinates, "\
                                      "so the PDB writer has terminated. "
                raise ValueError(print_error_message)

            PDB_atom_line_format = ('ATOM  %5s %-4s%1s%-4s%1s%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f      %-4s%2s%-2s\n')

            atom_alternate_location_List = []
            residue_code_insertion_List = []
            x_List = []
            y_List = []
            z_List = []
            atom_occupancy_List = []
            atom_bfactor_List = []
            Element_List = []

            # lock occupany factor at 1 (instead of: atom.occupancy)
            locked_occupany_factor = 1.00
            Max_No_atoms_in_base10 = 99999  # 99,999 for atoms in psf/pdb


            Res_No_iteration_corrected_List =[]
            Res_Chain_iteration_corrected_List = []
            residue_ID_list = []
            for i, atom_iter in enumerate( stuct_only_iteration.atoms):
                residue_ID_int = int(unique_residue_data_dict[residue_data_list[i]])
                Res_Chain_iteration_corrected_List.append(base10_to_base26_alph(int(residue_ID_int
                                                                                        / (self.Max_Residue_No + 1)))[-1:]
                                                          )
                Res_ID_adder = int((residue_ID_int % self.Max_Residue_No) % (self.Max_Residue_No))
                if int(Res_ID_adder) == 0:
                    Res_No_iteration_corrected_List.append(int(self.Max_Residue_No))
                else:
                    Res_No_iteration_corrected_List.append(Res_ID_adder)

                residue_ID_list.append(residue_ID_int)

            # This converts the atom name in the CHARMM psf and pdb files to unique atom names
            unique_Individual_atom_names_dict, \
            Individual_atom_names_List, \
            Missing_Bead_to_atom_name = unique_atom_naming(stuct_only_iteration, residue_ID_list, residue_names_list,
                                                           bead_to_atom_name_dict=self.bead_to_atom_name_dict)

            if None in [unique_Individual_atom_names_dict, Individual_atom_names_List, Missing_Bead_to_atom_name]:
                self.input_error = True
                print_error_message = 'ERROR: The unique_atom_naming function failed while '\
                                      'running the charmm_writer function. Ensure the proper inputs are ' \
                                      'in the bead_to_atom_name_dict.'

                raise ValueError(print_error_message)


            for  coord_iter, atom_coordinates in enumerate(All_atom_coordinates):

                for PDB_residue_count in  stuct_only_iteration.residues:
                    Segment_ID = ''
                    atom_iteration = sorted(PDB_residue_count.atoms, key=lambda atom: atom.number)
                    for atom_iteration_2 in atom_iteration:
                        x, y, z = atom_coordinates[atom_iteration_2.idx]
                        atom_alternate_location_List.append(atom_iteration_2.altloc)
                        residue_code_insertion_List.append(PDB_residue_count.insertion_code[:1])
                        x_List.append(x)
                        y_List.append(y)
                        z_List.append(z)
                        atom_occupancy_List.append(atom_iteration_2.occupancy)
                        atom_bfactor_List.append(atom_iteration_2.bfactor)
                        Element_List.append(Element[atom_iteration_2.atomic_number].upper())

                for v, atom_iter_1 in enumerate( stuct_only_iteration.atoms):

                    if v + 1 > Max_No_atoms_in_base10:
                        atom_number = base10_to_base16_alph_num(v + 1)

                    else:
                        atom_number = v + 1

                    output_write.write(PDB_atom_line_format % (atom_number, Individual_atom_names_List[v],
                                                               atom_alternate_location_List[v],
                                                               str(residue_names_list[v])[:self.No_1st_values_res_name],
                                                               Res_Chain_iteration_corrected_List[v],
                                                               Res_No_iteration_corrected_List[v],
                                                               residue_code_insertion_List[v],
                                                               x_List[v], y_List[v], z_List[v],
                                                               locked_occupany_factor, fix_atoms_list[v],
                                                               Segment_ID, Element_List[v], ''))

                output_write.write('%-80s\n' % 'END')

            output_write.close()

            # **********************************
            # **********************************
            # pdb writer (end)
            # **********************************
            # **********************************
