import os
from warnings import warn
from xml.dom import minidom

import parmed as pmd

from mbuild.utils.io import has_foyer

import mbuild as mb


def specific_ff_to_residue(structure,
                           forcefield_selection=None,
                           residues=None,
                           reorder_res_in_pdb_psf=False,
                           box=None,
                           boxes_for_simulation=1):
    """ToDo: List What this function does??
    Parameters
    ----------
    structure: mb.Compound
        The mBuild Compound object
    forcefield_selection: str or dictionary, default=None
        Apply a forcefield to the output file by selecting a force field xml file with
        its path or by using the standard force field name provided the `foyer` package.
        See Notes for further details on this argument
    residues: str of list of str, default=None
        Labels of unique residues in the Compound. Residues are assigned by
        checking against Compound.name.  Only supply residue names as 3 characters
        strings, as the residue names are truncated to 3 characters to fit in the
        psf and pdb file.
    reorder_res_in_pdb_psf: bool, default=False
        ToDo: Add Description fot this argument
    box: list, default=None
        list of 3 positive float values or the dimensions [x, y ,z]
        for structure in nanometers (nm). This is to add/override or change the structures dimensions.
        Ex: [1,2,3]
    boxes_for_simulation: {1, 2}
        Gibbs or grand canonical ensembles are examples of where the boxes_for_simulation would be 2

    Notes
    -----
    To write the NAMD/GOMC force field, pdb, and psf files, the
    residues and forcefields must be provided in a str or
    dictionary.  If a dictionary is provided all residues must
    be specified to a force field.

    Example dict for FF file: {'ETH' : 'oplsaa.xml', 'OCT': 'path_to file/trappe-ua.xml'}
    Example str for FF file: 'path_to file/trappe-ua.xml'
    Example dict for standard FF names : {'ETH' : 'oplsaa', 'OCT': 'trappe-ua'}
    Example str for standard FF names: 'trappe-ua'
    Example of a mixed dict with both : {'ETH' : 'oplsaa', 'OCT': 'path_to file/'trappe-ua.xml'}

    Returns
    -------
    structure: parmed.Structure
        parmed structure with applied force field
    coulomb14scaler_dict: dict
        a dictionary with the 1,4-colombic scalers for each residue
            (i.e., a different force field could on each residue)
    lj14_scaler_dict: dict
        a dictionary with the 1,4-LJ scalers for each residue
        (i.e., a different force field could on each residue)
    residues_applied_list: list
        list of residues (i.e., list of stings).
        These are all the residues in which the force field actually applied
    """

    if has_foyer:
        from foyer import Forcefield
        from foyer.forcefields import forcefields
    else:
        warn(
            'Package foyer is not installed. '
            'Please install it using conda install -c conda-forge foyer'
        )
        return None, None, None, None

    if forcefield_selection is None:
        warn('Please enter either the forcefields for the forcefield_selection variable')
        return None, None, None, None

    elif forcefield_selection is not None and not isinstance(forcefield_selection, dict):
        warn(
            'The force field selection (forcefield_selection) '
            'is not a dictionary. Please enter a dictionary '
            'with all the residues specified to a force field '
            '-> Ex: {"Water" : "oplsaa", "OCT": "path/trappe-ua.xml"}, '
            'Note: the file path must be specified the force field file '
            'or by using the standard force field name provided the `foyer` package.'
        )
        return None, None, None, None

    if residues is None or isinstance(residues, list) == False:
        print('Please enter the residues in the Specific_FF_to_residue function')
        return None, None, None, None


    if not isinstance(reorder_res_in_pdb_psf, bool):
        print(
            'Please enter the reorder_res_in_pdb_psf '
            'in the Specific_FF_to_residue function (i.e., True or False)'
        )
        return None, None, None, None

    if box is not None:
        box_ang = []
        box_length = len(box)
        if box_length != 3:
            warn('Please enter all 3 values for the box dimensions.')
            return None, None, None, None
        for box_iter in range(0, len(box)):
            if isinstance(box[box_iter], str):
                warn('Please enter all positive or 0 values for the box dimensions.')
                return None, None, None, None
            if box[box_iter] < 0:
                warn('Please enter all positive or 0 values for the box dimensions.')
                return None, None, None, None
            # change from nm to Angstroms
            box_ang.append(box[box_iter] * 10)

    if boxes_for_simulation not in [1, 2]:
        warn('Please enter boxes_for_simulation equal the integer 1 or 2.')
        return None, None, None, None



    forcefield_keys_list = []
    if forcefield_selection is not None:
        for res in forcefield_selection.keys():
            forcefield_keys_list.append(res)
        ff_data = forcefield_selection


    if forcefield_keys_list == [] and len(residues) != 0 :
        warn(
            'The forcefield_selection variable are not provided, but there are residues provided.'
        )
        return None, None, None, None
    elif forcefield_keys_list != [] and len(residues) == 0 :
        warn(
            'The residues variable is and empty list but there are forcefield_selection variables provided.'
        )
        return None, None, None, None


    user_entered_ff_with_path_dict = {}  # True means user entered the path, False is a standard foyer FF with no path
    for z in range(0, len(forcefield_keys_list)):
        for res_i in range(0, len(residues)):
            if residues[res_i] == forcefield_keys_list[z]:
                if os.path.splitext(ff_data[forcefield_keys_list[z]])[1] == '.xml' and len(residues) != 0:
                    user_entered_ff_with_path_dict.update({residues[res_i]: True})
                elif os.path.splitext(ff_data[forcefield_keys_list[z]])[1] == '' and len(residues) != 0:
                    user_entered_ff_with_path_dict.update({residues[res_i]: False})
                else:
                    warn(
                        'Please make sure are entering the correct '
                        'foyer FF name and not a path to a FF file. ' 
                        'If you are entering a path to a FF file, '
                        'please use the forcefield_files variable.'
                    )
                    return None, None, None, None


    coulomb14scaler_dict = {}
    lj14_scaler_dict = {}
    for j in range(0, len(forcefield_keys_list)):
        residue_iteration = forcefield_keys_list[j]
        if user_entered_ff_with_path_dict[residue_iteration] :
            ff_for_residue_iteration = ff_data[residue_iteration]

            try:
                read_xlm_iteration = minidom.parse(ff_for_residue_iteration)
            except:
                warn('Please make sure are entering the correct foyer FF path, including the FF file name.xml ' +
                     'If you are using the pre-build FF files in foyer, please us the forcefield_names variable.')
                return None, None, None, None

        elif not user_entered_ff_with_path_dict[residue_iteration]:
            ff_for_residue_iteration = ff_data[residue_iteration]
            ff_names_path_iteration = forcefields.get_ff_path()[0] + '/xml/' + ff_for_residue_iteration + '.xml'
            try:
                read_xlm_iteration = minidom.parse(ff_names_path_iteration)
            except:
                warn('Please make sure are entering the correct foyer FF name and not a path to a FF file.' +
                     'If you are entering a path to a FF file, please us the forcefield_files variable')
                return None, None, None, None
        lj_coul_1_4_values = read_xlm_iteration.getElementsByTagName("NonbondedForce")

        for Scaler in lj_coul_1_4_values:
            coulomb14scaler_dict.update({residue_iteration: float(Scaler.getAttribute("coulomb14scale"))})
            lj14_scaler_dict.update({residue_iteration: float(Scaler.getAttribute("lj14scale"))})

    # calculate the initial number of atoms for later comparison
    initial_no_atoms = len(structure.to_parmed().atoms)

    # Check to see if it is an empty mbuild.Compound and set intial atoms to 0
    # note empty mbuild.Compound will read 1 atoms but there is really noting there
    if str(structure.to_parmed()) == '<Structure 1 atoms; 1 residues; 0 bonds; PBC (orthogonal); NOT parametrized>' \
           and str(structure.children) == 'OrderedSet()' and initial_no_atoms == 1 \
           and str(structure.pos)=='[0. 0. 0.]' \
           and str(structure.to_parmed().atoms[0]) == '<Atom Compound [0]; In RES 0>' \
           and str(structure.to_parmed().residues[0]) ==  '<Residue RES[0]>' \
            and len(structure.children) == 0 :
        # there are no real atoms so set initial atoms to 0
        initial_no_atoms = 0


    # add the FF to the residues
    compound_box_infor = structure.to_parmed(residues=residues)
    new_structure = pmd.Structure()
    new_structure.box = compound_box_infor.box

    # prepare all compound and remove nested compounds
    no_layers_to_check_for_residues = 4

    for j in range(0, no_layers_to_check_for_residues):
        new_compound_iter = mb.Compound()
        new_compound_iter.periodicity[0] = structure.periodicity[0]
        new_compound_iter.periodicity[1] = structure.periodicity[1]
        new_compound_iter.periodicity[2] = structure.periodicity[2]
        if structure.name in residues:
            if len(structure.children) == 0:
                warn('Warning: This residue is the atom, and is a single atom., ' + str(structure.name))
                new_compound_iter.add(mb.compound.clone(structure))

            elif len(structure.children) > 0:

                new_compound_iter.add(mb.compound.clone(structure))

        else:
            for child in structure.children:
                if len(child.children) == 0:
                    if child.name not in residues:
                        warn('ERROR: All the residues are not specified')
                        return None, None, None, None

                    else:
                        new_compound_iter.add(mb.compound.clone(child))

                elif len(child.children) > 0:
                    if child.name in residues:
                        new_compound_iter.add(mb.compound.clone(child))
                    else:
                        for sub_child in child.children:
                            if sub_child.name in residues:
                                new_compound_iter.add(mb.compound.clone(sub_child))
                            else:
                                if len(sub_child.children) == 0 and (child.name not in residues):
                                    warn('ERROR: All the residues are not specified')
                                    return None, None, None, None

        structure = new_compound_iter

    residues_applied_list = []
    residue_orig_order_list = []
    for child in structure.children:
        if child.name not in residue_orig_order_list:
            residue_orig_order_list.append(child.name)
    for res_reorder_iter in range(0, len(residues)):
        if residues[res_reorder_iter] not in residue_orig_order_list:
            text_to_print_1 = "All the residues were not used from the forcefield_selection " + \
                              "string or dictionary.  There may be residues below other specified residues " + \
                              "in the mbuild.Compound hierarchy.  If so, all the highest listed residues pass " + \
                              "down the force fields through the hierarchy.  Alternatively, " + \
                              "residues that are not in the structure may have been specified. "
            text_to_print_2 = "Note: This warning will appear if you are using the CHARMM pdb and psf writers " + \
                              "2 boxes, and the boxes do not contain all the residues in each box."
            if boxes_for_simulation == 1:
                warn(text_to_print_1)
                return None, None, None, None
            if boxes_for_simulation == 2:
                warn(text_to_print_1 + text_to_print_2)

    if not reorder_res_in_pdb_psf:
        residues = residue_orig_order_list
    elif reorder_res_in_pdb_psf:
        print("Information: the output file are being reordered in via the residues list's sequence. ")

    for i in range(0, len(residues)):
        children_in_iteration = False
        new_compound_iteration = mb.Compound()
        new_compound_iter.periodicity[0] = structure.periodicity[0]
        new_compound_iter.periodicity[1] = structure.periodicity[1]
        new_compound_iter.periodicity[2] = structure.periodicity[2]
        new_structure_iteration = pmd.Structure()
        new_structure_iteration.box = compound_box_infor.box
        for child in structure.children:
            if ff_data.get(child.name) is None:
                warn('ERROR: All residues are not specified in the force_field dictionary')
                return None, None, None, None

            if child.name == residues[i]:
                children_in_iteration = True
                new_compound_iteration.add(mb.compound.clone(child))

        if children_in_iteration:
            if user_entered_ff_with_path_dict[residues[i]]:
                ff_iteration = Forcefield(ff_data[residues[i]])
                residues_applied_list.append(residues[i])
            elif not user_entered_ff_with_path_dict[residues[i]]:
                ff_iteration = Forcefield(name=ff_data[residues[i]])
                residues_applied_list.append(residues[i])

            new_structure_iteration = ff_iteration.apply(new_compound_iteration, residues=[residues[i]])
            new_structure = new_structure + new_structure_iteration

    if box is not None:
        new_structure.box[0] = box_ang[0]
        new_structure.box[1] = box_ang[1]
        new_structure.box[2] = box_ang[2]

    structure = new_structure

    # calculate the final number of atoms
    final_no_atoms = len(structure.atoms)

    if final_no_atoms != initial_no_atoms:
        warn('ERROR: The initial number of atoms send to the force field analysis is not the same ' +
             'as the final number of atoms analyzed.  The intial number of atoms was ' + str(initial_no_atoms) +
             ' and the final number of atoms was ' + str(final_no_atoms) + '. Please ensure that all ' +
             'the residues names that are in the initial Compound are listed in the ' +
             'residues list (i.e., the residues variable)')
        return None, None, None, None

    return structure, coulomb14scaler_dict, lj14_scaler_dict, residues_applied_list
