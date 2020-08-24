from foyer import Forcefield
from foyer.forcefields import forcefields
import parmed as pmd
import mbuild as mb
import xml.dom.minidom
from warnings import warn

def Specific_FF_to_residue(structure , forcefield_files= None, forcefield_names= None,
                           residues= None, reorder_res_in_pdb_psf= False):

    #input:
        # structure =  compound structure
        #  forcefield_files =  dictionary of residues to force fields.
            # Ex: {'Water' : 'oplsaa.xml', 'OCT': 'trappe-ua.xml'}
        #  forcefield_names =  dictionary of residues to force fields.
            # Ex: {'Water' : 'oplsaa', 'OCT': 'trappe-ua'}
        # Note: either the forcefield_files or forcefield_names can be provided, not both.
                # the forcefield_files or forcefield_names can be not be mixed combination.
        # residues = list of residues
        # note FFs and residues must be in sequenctial order

    #Returns:
        # structure = parmed with applied force field
        # coulomb14scaler_dict = a dictionary with the 1,4-colombic scalers for each residue
                                # (i.e., a different force field could on each residue)
        # LJ14scaler_dict  = a dictionary with the 1,4-LJ scalers for each residue
                            # (i.e., a different force field could on each residue)

    if forcefield_names is None and forcefield_files is None:
        return warn('Please enter either the forcefield_files or forcefield_names, neither were provided')

    if forcefield_names != None and forcefield_files != None:
        return warn('Please enter either the forcefield_files or forcefield_names, not both')

    elif forcefield_files != None and forcefield_names is None and not isinstance(forcefield_files, dict):
        return warn('The force field file (forcefield_files) is not a dictionary. Please enter a dictionary'+
                      "with all the residues specified to a force field"+
                     "-> Ex: {'Water' : 'oplsaa.xml', 'OCT': 'trappe-ua.xml'}, "+
                     "Note: the file path must be specified the force field file")

    elif forcefield_names != None and forcefield_files is None and not isinstance(forcefield_names, dict) and not isinstance(forcefield_names, str):
        return warn('The force field names (forcefield_names) is not a string or a dictionary with' +
                    ' all the residues specified to a force field.' +
                    "-> String Ex: 'trappe-ua' or 'oplsaa'  ." +
                    "Otherise provided a dictionary with all the residues specified to a force field " +
                    "->Dictionary Ex: {'Water' : 'oplsaa', 'OCT': 'trappe-ua'}, " +
                    "Note: the file path must be specified the force field file")


    if residues is None:
        print('please enter the residues in the Specific_FF_to_residue function')
    if reorder_res_in_pdb_psf is None:
        print('please enter the reorder_res_in_pdb_psf in the Specific_FF_to_residue function')


    forcefield_keys_list = []
    Use_FF_files = False
    Use_FF_names = False
    if forcefield_names != None and forcefield_files is None:
        Use_FF_names = True
        for res in forcefield_names.keys():
            forcefield_keys_list.append(res)
        FF_data = forcefield_names
    elif forcefield_files != None and forcefield_names is None:
        Use_FF_files = True
        for res in forcefield_files.keys():
            forcefield_keys_list.append(res)
        FF_data = forcefield_files

    coulomb14scaler_dict = {}
    LJ14scaler_dict = {}
    for j in range(0, len(forcefield_keys_list)):
        if Use_FF_files == True:
            residue_iteration = forcefield_keys_list[j]
            FF_for_residue_iteration = FF_data[residue_iteration]
            try:
                read_xlm_iteration = xml.dom.minidom.parse(FF_for_residue_iteration)
            except:
                return warn('Please make sure are entering the correct foyer FF path, including the FF file name.xml ' +
                            'If you are using the pre-build FF files in foyer, please us the forcefield_names variable.')
            LJ_Coul_1_4_values = read_xlm_iteration.getElementsByTagName("NonbondedForce")

        elif Use_FF_names == True:
            residue_iteration = forcefield_keys_list[j]
            FF_for_residue_iteration = FF_data[residue_iteration]
            FF_names_path_iteration = forcefields.get_ff_path()[0] +'/xml/'+FF_for_residue_iteration+'.xml'
            try:
                read_xlm_iteration = xml.dom.minidom.parse(FF_names_path_iteration)
            except:
                return warn('Please make sure are entering the correct foyer FF name and not a path to a FF file.' +
                            'If you are entering a path to a FF file, please us the forcefield_files variable')
            LJ_Coul_1_4_values = read_xlm_iteration.getElementsByTagName("NonbondedForce")

        for Scaler in LJ_Coul_1_4_values:
            coulomb14scaler_dict.update({residue_iteration: float(Scaler.getAttribute("coulomb14scale"))})
            LJ14scaler_dict.update({residue_iteration: float(Scaler.getAttribute("lj14scale"))})





    # add the FF to the residues
    compound_box_infor = structure.to_parmed(residues=residues)
    new_structure = pmd.Structure()
    new_structure.box = compound_box_infor.box



    #prepare all compound and remove nested compounds
    restructure_compound = False
    for child in structure.children:
        if child.name not in residues:
            restructure_compound = True

    No_layers_to_check_for_residues = 4
    if restructure_compound == True:
        for j in range(0, No_layers_to_check_for_residues):
            new_compound_iter = mb.Compound()
            new_compound_iter.periodicity[0] = structure.periodicity[0]
            new_compound_iter.periodicity[1] = structure.periodicity[1]
            new_compound_iter.periodicity[2] = structure.periodicity[2]

            for child in structure.children:
                if len(child.children)==0:
                    if child.name not in residues:
                        print('All the residues are not specified')
                        return warn('All the residues are not specified')

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
                                if len(sub_child.children)==0 and (child.name not in residues):
                                    print('All the residues are not specified')
                                    return warn('All the residues are not specified')

            structure = new_compound_iter


    Residue_orig_order_list = []
    for child in structure.children:
        if child.name not in  Residue_orig_order_list:
            Residue_orig_order_list.append(child.name)


    if reorder_res_in_pdb_psf==False:
        residues = Residue_orig_order_list
    elif reorder_res_in_pdb_psf==True:
        residues = residues
    else:
        print("ERROR residues = Residue_orig_order_list or residues= residues not properly specified ")
        return warn("ERROR residues = Residue_orig_order_list or residues= residues not properly specified ")

    for i in range(0, len(residues)):
        children_in_iteration = False
        new_compound_iteration = mb.Compound()
        new_compound_iteration.periodicity[0] = structure.periodicity[0]
        new_compound_iteration.periodicity[1] = structure.periodicity[1]
        new_compound_iteration.periodicity[2] = structure.periodicity[2]
        new_structure_iteration = pmd.Structure()
        new_structure_iteration.box = compound_box_infor.box
        for child in structure.children:
            if FF_data.get(child.name) is None:
                print('All residues are not specified in the force_field dictionary')
                return warn('All residues are not specified in the force_field dictionary')

            if child.name == residues[i]:
                children_in_iteration = True
                new_compound_iteration.add(mb.compound.clone(child))

        if children_in_iteration == True:
            if Use_FF_files == True:
                FF_iteration = Forcefield(FF_data[residues[i]])
            elif Use_FF_names == True:
                FF_iteration = Forcefield(name=FF_data[residues[i]])

            new_structure_iteration = FF_iteration.apply(new_compound_iteration, residues=[residues[i]])
            new_structure = new_structure + new_structure_iteration

    structure = new_structure
    return structure, coulomb14scaler_dict, LJ14scaler_dict

