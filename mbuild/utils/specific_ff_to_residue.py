import os
from warnings import warn
from xml.dom import minidom

import parmed as pmd

import mbuild as mb
from mbuild.compound import Compound
from mbuild.utils.io import has_foyer


def specific_ff_to_residue(
    structure,
    forcefield_selection=None,
    residues=None,
    reorder_res_in_pdb_psf=False,
    boxes_for_simulation=1,
):

    """
    Takes the mbuild Compound or mbuild Box structure and applies the selected
    force field to the corresponding residue via foyer.
    Note: a residue is defined as a molecule in this case, so it is not
    designed for applying a force field to a protein.

    Parameters
    ----------
    structure: mbuild Compound object or mbuild Box object;
        The mBuild Compound object or mbuild Box object, which contains the molecules
        (or empty box) that will have the force field applied to them.
    forcefield_selection: str or dictionary, default=None
        Apply a forcefield to the output file by selecting a force field xml file with
        its path or by using the standard force field name provided the `foyer` package.
        Example dict for FF file: {'ETH' : 'oplsaa.xml', 'OCT': 'path_to file/trappe-ua.xml'}
        Example str for FF file: 'path_to file/trappe-ua.xml'
        Example dict for standard FF names : {'ETH' : 'oplsaa', 'OCT': 'trappe-ua'}
        Example str for standard FF names: 'trappe-ua'
        Example of a mixed dict with both : {'ETH' : 'oplsaa', 'OCT': 'path_to file/'trappe-ua.xml'}
    residues: list, [str, ..., str], default=None
        Labels of unique residues in the Compound. Residues are assigned by
        checking against Compound.name.  Only supply residue names as 4 characters
        strings, as the residue names are truncated to 4 characters to fit in the
        psf and pdb file.
    reorder_res_in_pdb_psf: bool, default=False
        This option provides the ability to reorder the residues/molecules from the original
        structure's order.  If True, the residues will be reordered as they appear in the residues
        variable.  If False, the order will be the same as entered in the original structure.
    boxes_for_simulation: int [1, 2], default = 1
        Gibbs (GEMC) or grand canonical (GCMC) ensembles are examples of where the boxes_for_simulation would be 2.
        Canonical (NVT) or isothermal–isobaric (NPT) ensembles are example with the boxes_for_simulation equal to 1.
        Note: the only valid options are 1 or 2.

    Returns
    -------
    list, [structure, coulomb14scalar_dict, lj14_scalar_dict, residues_applied_list]
        structure: parmed.Structure
            parmed structure with applied force field
        coulomb14scalar_dict: dict
            a dictionary with the 1,4-colombic scalars for each residue
                (i.e., a different force field could on each residue)
        lj14_scalar_dict: dict
            a dictionary with the 1,4-LJ scalars for each residue
            (i.e., a different force field could on each residue)
        residues_applied_list: list
            list of residues (i.e., list of stings).
            These are all the residues in which the force field actually applied

    Notes
    -----
    To write the NAMD/GOMC force field, pdb, psf, and force field
    (.inp) files, the residues and forcefields must be provided in
    a str or dictionary. If a dictionary is provided all residues must
    be specified to a force field if the boxes_for_simulation is equal to 1.

    Generating an empty box (i.e., pdb and psf files):
    Enter residues = [], but the accompanying structure must be an empty mb.Box.
    However, when doing this, the forcefield_selection must be supplied,
    or it will provide an error (i.e., forcefield_selection can not be equal to None).

    In this current FF/psf/pdb writer, a residue type is essentially a molecule type.
    Therefore, it can only correctly write systems where every bead/atom in the molecule
    has the same residue name, and the residue name is specific to that molecule type.
    For example: a protein molecule with many residue names is not currently supported,
    but is planned to be supported in the future.
    """

    if has_foyer:
        from foyer import Forcefield
        from foyer.forcefields import forcefields
    else:
        print_error_message = (
            "Package foyer is not installed. "
            "Please install it using conda install -c conda-forge foyer"
        )
        raise ImportError(print_error_message)

    if not isinstance(structure, (Compound, mb.Box)):
        print_error_message = (
            "ERROR: The structure expected to be of type: "
            "{} or {}, received: {}".format(
                type(Compound()),
                type(mb.Box(lengths=[1, 1, 1])),
                type(structure),
            )
        )
        raise TypeError(print_error_message)

    print("forcefield_selection = " + str(forcefield_selection))
    if forcefield_selection is None:
        print_error_message = (
            "Please the force field selection (forcefield_selection) as a dictionary "
            "with all the residues specified to a force field "
            '-> Ex: {"Water" : "oplsaa", "OCT": "path/trappe-ua.xml"}, '
            "Note: the file path must be specified the force field file "
            "or by using the standard force field name provided the `foyer` package."
        )
        raise TypeError(print_error_message)

    elif forcefield_selection is not None and not isinstance(
        forcefield_selection, dict
    ):
        print_error_message = (
            "The force field selection (forcefield_selection) "
            "is not a dictionary. Please enter a dictionary "
            "with all the residues specified to a force field "
            '-> Ex: {"Water" : "oplsaa", "OCT": "path/trappe-ua.xml"}, '
            "Note: the file path must be specified the force field file "
            "or by using the standard force field name provided the `foyer` package."
        )
        raise TypeError(print_error_message)

    if residues is None or not isinstance(residues, list):
        print_error_message = (
            "Please enter the residues in the Specific_FF_to_residue function."
        )
        raise TypeError(print_error_message)

    if not isinstance(reorder_res_in_pdb_psf, bool):
        print_error_message = (
            "Please enter the reorder_res_in_pdb_psf "
            "in the Specific_FF_to_residue function (i.e., True or False)."
        )
        raise TypeError(print_error_message)

    print_error_message_for_boxes_for_simulatiion = (
        "ERROR: Please enter boxes_for_simulation equal " "the integer 1 or 2."
    )
    if not isinstance(boxes_for_simulation, int):
        raise TypeError(print_error_message_for_boxes_for_simulatiion)

    elif isinstance(boxes_for_simulation, int) and boxes_for_simulation not in [
        1,
        2,
    ]:
        raise ValueError(print_error_message_for_boxes_for_simulatiion)

    forcefield_keys_list = []
    if forcefield_selection is not None:
        for res in forcefield_selection.keys():
            forcefield_keys_list.append(res)
        ff_data = forcefield_selection

    if forcefield_keys_list == [] and len(residues) != 0:
        print_error_message = "The forcefield_selection variable are not provided, but there are residues provided."
        raise ValueError(print_error_message)

    elif forcefield_keys_list != [] and len(residues) == 0:
        print_error_message = (
            "The residues variable is an empty list but there are "
            "forcefield_selection variables provided."
        )
        raise ValueError(print_error_message)

    user_entered_ff_with_path_dict = (
        {}
    )  # True means user entered the path, False is a standard foyer FF with no path
    for z in range(0, len(forcefield_keys_list)):
        for res_i in range(0, len(residues)):
            if residues[res_i] == forcefield_keys_list[z]:
                if (
                    os.path.splitext(ff_data[forcefield_keys_list[z]])[1]
                    == ".xml"
                    and len(residues) != 0
                ):
                    user_entered_ff_with_path_dict.update(
                        {residues[res_i]: True}
                    )
                elif (
                    os.path.splitext(ff_data[forcefield_keys_list[z]])[1] == ""
                    and len(residues) != 0
                ):
                    user_entered_ff_with_path_dict.update(
                        {residues[res_i]: False}
                    )
                else:
                    print_error_message = (
                        r"Please make sure you are entering the correct "
                        "foyer FF name and not a path to a FF file. "
                        "If you are entering a path to a FF file, "
                        "please use the forcefield_files variable with the "
                        "proper XML extension (.xml)."
                    )
                    raise ValueError(print_error_message)

    coulomb14scalar_dict = {}
    lj14_scalar_dict = {}
    for j in range(0, len(forcefield_keys_list)):
        residue_iteration = forcefield_keys_list[j]
        if user_entered_ff_with_path_dict[residue_iteration]:
            ff_for_residue_iteration = ff_data[residue_iteration]
            try:
                read_xlm_iteration = minidom.parse(ff_for_residue_iteration)

            except:
                print_error_message = (
                    "Please make sure you are entering the correct foyer FF path, "
                    "including the FF file name.xml "
                    "If you are using the pre-build FF files in foyer, "
                    "only use the string name without any extension."
                )
                raise ValueError(print_error_message)
        elif not user_entered_ff_with_path_dict[residue_iteration]:
            ff_for_residue_iteration = ff_data[residue_iteration]
            ff_names_path_iteration = (
                forcefields.get_ff_path()[0]
                + "/xml/"
                + ff_for_residue_iteration
                + ".xml"
            )
            try:
                read_xlm_iteration = minidom.parse(ff_names_path_iteration)
            except:
                print_error_message = (
                    "Please make sure you are entering the correct foyer FF name, or the "
                    "correct file extension (i.e., .xml, if required)."
                )
                raise ValueError(print_error_message)
        lj_coul_1_4_values = read_xlm_iteration.getElementsByTagName(
            "NonbondedForce"
        )

        for Scalar in lj_coul_1_4_values:
            coulomb14scalar_dict.update(
                {
                    residue_iteration: float(
                        Scalar.getAttribute("coulomb14scale")
                    )
                }
            )
            lj14_scalar_dict.update(
                {residue_iteration: float(Scalar.getAttribute("lj14scale"))}
            )

    # Check to see if it is an empty mbuild.Compound and set intial atoms to 0
    # note empty mbuild.Compound will read 1 atoms but there is really noting there
    if isinstance(structure, Compound):
        if len(structure.children) == 0:
            # there are no real atoms in the Compound so the test fails. User should use mbuild.Box
            print_error_message = (
                "ERROR: If you are not providing an empty box, "
                "you need to specify the atoms/beads as children in the mb.Compound. "
                "If you are providing and empty box, please do so by specifying and "
                "mbuild Box ({})".format(type(mb.Box(lengths=[1, 1, 1])))
            )
            raise TypeError(print_error_message)
        else:
            initial_no_atoms = len(structure.to_parmed().atoms)

    # calculate the initial number of atoms for later comparison
    if isinstance(structure, mb.Box):
        lengths = structure.lengths
        angles = structure.angles

        structure = mb.Compound()
        structure.box = mb.Box(lengths=lengths, angles=angles)
        initial_no_atoms = 0

    # add the FF to the residues
    compound_box_infor = structure.to_parmed(residues=residues)
    new_structure = pmd.Structure()
    new_structure.box = compound_box_infor.box

    # prepare all compound and remove nested compounds
    no_layers_to_check_for_residues = 3

    print_error_message_all_res_not_specified = (
        "ERROR: All the residues are not specified, or "
        "the residues entered does not match the residues that "
        "were found and built for structure."
    )
    for j in range(0, no_layers_to_check_for_residues):
        new_compound_iter = mb.Compound()
        new_compound_iter.periodicity = structure.periodicity
        if structure.name in residues:
            if len(structure.children) == 0:
                warn(
                    "Warning: This residue is the atom, and is a single atom., "
                    + str(structure.name)
                )
                new_compound_iter.add(mb.compound.clone(structure))

            elif len(structure.children) > 0:

                new_compound_iter.add(mb.compound.clone(structure))

        else:
            for child in structure.children:
                if len(child.children) == 0:
                    if child.name not in residues:
                        raise ValueError(
                            print_error_message_all_res_not_specified
                        )

                    else:
                        new_compound_iter.add(mb.compound.clone(child))

                elif len(child.children) > 0:
                    if child.name in residues:
                        new_compound_iter.add(mb.compound.clone(child))
                    else:
                        for sub_child in child.children:
                            if sub_child.name in residues:
                                new_compound_iter.add(
                                    mb.compound.clone(sub_child)
                                )

                            else:
                                if len(sub_child.children) == 0 and (
                                    child.name not in residues
                                ):

                                    raise ValueError(
                                        print_error_message_all_res_not_specified
                                    )

        structure = new_compound_iter

    residues_applied_list = []
    residue_orig_order_list = []
    for child in structure.children:
        if child.name not in residue_orig_order_list:
            residue_orig_order_list.append(child.name)
    for res_reorder_iter in range(0, len(residues)):
        if residues[res_reorder_iter] not in residue_orig_order_list:
            text_to_print_1 = (
                "All the residues were not used from the forcefield_selection "
                "string or dictionary. There may be residues below other "
                "specified residues in the mbuild.Compound hierarchy. "
                "If so, all the highest listed residues pass down the force "
                "fields through the hierarchy. Alternatively, residues that "
                "are not in the structure may have been specified. "
            )
            text_to_print_2 = (
                "Note: This warning will appear if you are using the CHARMM pdb and psf writers "
                + "2 boxes, and the boxes do not contain all the residues in each box."
            )
            if boxes_for_simulation == 1:
                warn(text_to_print_1)
                raise ValueError(text_to_print_1)
            if boxes_for_simulation == 2:
                warn(text_to_print_1 + text_to_print_2)

    if not reorder_res_in_pdb_psf:
        residues = residue_orig_order_list
    elif reorder_res_in_pdb_psf:
        print(
            "INFO: the output file are being reordered in via the residues list's sequence."
        )

    for i in range(0, len(residues)):
        children_in_iteration = False
        new_compound_iteration = mb.Compound()
        new_compound_iter.periodicity = structure.periodicity
        new_structure_iteration = pmd.Structure()
        new_structure_iteration.box = compound_box_infor.box
        for child in structure.children:
            if ff_data.get(child.name) is None:
                print_error_message = "ERROR: All residues are not specified in the force_field dictionary"
                raise ValueError(print_error_message)

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

            new_compound_iteration.box = None
            new_structure_iteration = ff_iteration.apply(
                new_compound_iteration, residues=[residues[i]]
            )
            new_structure = new_structure + new_structure_iteration

    structure = new_structure

    # calculate the final number of atoms
    final_no_atoms = len(structure.atoms)

    if final_no_atoms != initial_no_atoms:
        print_error_message = (
            "ERROR: The initial number of atoms sent to the force field analysis is "
            "not the same as the final number of atoms analyzed. "
            "The initial number of atoms was {} and the final number of atoms was {}. "
            "Please ensure that all the residues names that are in the initial "
            "Compound are listed in the residues list "
            "(i.e., the residues variable).".format(
                initial_no_atoms, final_no_atoms
            )
        )
        raise ValueError(print_error_message)

    return [
        structure,
        coulomb14scalar_dict,
        lj14_scalar_dict,
        residues_applied_list,
    ]
