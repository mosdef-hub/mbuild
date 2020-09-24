from foyer import Forcefield
from foyer.forcefields import forcefields
import parmed as pmd
import mbuild as mb
import xml.dom.minidom
from warnings import warn

def Specific_FF_to_residue(structure , forcefield_files= None, forcefield_names= None,
                           residues= None, reorder_res_in_pdb_psf= False, box = None,
                           boxes_for_simulation = 1):

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
        # box ; list of 3 positive float values or the dimensions [x, y ,z]
             #for structure_1 in nanometers (nm)
             #This is to add/override or change the structures dimenstions. Ex: [1,2,3]
        # boxes_for_simulation; int of 1 or 2, default =1.  Gibbs or grand canonical ensembles
            #are examples of where the boxes_for_simulation would be 2

    #Returns:
        # structure = parmed with applied force field
        # coulomb14scaler_dict = a dictionary with the 1,4-colombic scalers for each residue
                                # (i.e., a different force field could on each residue)
        # LJ14scaler_dict  = a dictionary with the 1,4-LJ scalers for each residue
                            # (i.e., a different force field could on each residue)
        # residues_applied_list =  list of residues (i.e., list of stings)
                            # these are all the residues in which the force field actually applied

    if forcefield_names is None and forcefield_files is None:
        warn('Please enter either the forcefield_files or forcefield_names, neither were provided')
        return None, None, None, None

    if forcefield_names != None and forcefield_files != None:
        warn('Please enter either the forcefield_files or forcefield_names, not both')
        return None, None, None, None

    elif forcefield_files != None and forcefield_names is None and not isinstance(forcefield_files, dict):
        warn('The force field file (forcefield_files) is not a dictionary. Please enter a dictionary'+
                      "with all the residues specified to a force field"+
                     "-> Ex: {'Water' : 'oplsaa.xml', 'OCT': 'trappe-ua.xml'}, "+
                     "Note: the file path must be specified the force field file")
        return None, None, None, None

    elif forcefield_names != None and forcefield_files is None and not isinstance(forcefield_names, dict) and not isinstance(forcefield_names, str):
        warn('The force field names (forcefield_names) is not a string or a dictionary with' +
                    ' all the residues specified to a force field.' +
                    "-> String Ex: 'trappe-ua' or 'oplsaa'  ." +
                    "Otherise provided a dictionary with all the residues specified to a force field " +
                    "->Dictionary Ex: {'Water' : 'oplsaa', 'OCT': 'trappe-ua'}, " +
                    "Note: the file path must be specified the force field file")
        return None, None, None, None


    if residues is None:
        print('Please enter the residues in the Specific_FF_to_residue function')
    if reorder_res_in_pdb_psf is None:
        print('Please enter the reorder_res_in_pdb_psf in the Specific_FF_to_residue function')

    if box !=None :
        box_Ang = []
        box_length = len(box)
        if box_length != 3:
            warn('Please enter all 3 values for the box dimensions.')
            return None, None, None, None
        for box_iter in range(0, len(box)):
            if isinstance(box[ box_iter], str)==True:
                warn('Please enter all positive or 0 values for the box dimensions.')
                return None, None, None, None
            if box[ box_iter] < 0:
                warn('Please enter all positive or 0 values for the box dimensions.')
                return None, None, None, None
            # change from nm to Angstroms
            box_Ang.append(box[ box_iter]*10)

    if boxes_for_simulation not in [1, 2]:
        warn('Please enter boxes_for_simulation equal the interger 1 or 2.')
        return None, None, None, None

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
                warn('Please make sure are entering the correct foyer FF path, including the FF file name.xml ' +
                        'If you are using the pre-build FF files in foyer, please us the forcefield_names variable.')
                return None, None, None, None

        elif Use_FF_names == True:
            residue_iteration = forcefield_keys_list[j]
            FF_for_residue_iteration = FF_data[residue_iteration]
            FF_names_path_iteration = forcefields.get_ff_path()[0] +'/xml/'+FF_for_residue_iteration+'.xml'
            try:
                read_xlm_iteration = xml.dom.minidom.parse(FF_names_path_iteration)
            except:
                warn('Please make sure are entering the correct foyer FF name and not a path to a FF file.' +
                        'If you are entering a path to a FF file, please us the forcefield_files variable')
                return None, None, None, None
        LJ_Coul_1_4_values = read_xlm_iteration.getElementsByTagName("NonbondedForce")

        for Scaler in LJ_Coul_1_4_values:
            coulomb14scaler_dict.update({residue_iteration: float(Scaler.getAttribute("coulomb14scale"))})
            LJ14scaler_dict.update({residue_iteration: float(Scaler.getAttribute("lj14scale"))})





    # add the FF to the residues
    compound_box_infor = structure.to_parmed(residues=residues)
    new_structure = pmd.Structure()
    new_structure.box = compound_box_infor.box


    #prepare all compound and remove nested compounds
    No_layers_to_check_for_residues = 4

    for j in range(0, No_layers_to_check_for_residues):
        new_compound_iter = mb.Compound()
        new_compound_iter.periodicity[0] = structure.periodicity[0]
        new_compound_iter.periodicity[1] = structure.periodicity[1]
        new_compound_iter.periodicity[2] = structure.periodicity[2]
        if structure.name in residues:
            if len(structure.children) == 0:
                warn('ERROR: There are no atoms in this residue, '+str(structure.name))
                return None, None, None, None

            elif len(structure.children) > 0:

                new_compound_iter.add(mb.compound.clone(structure))

        else:
            for child in structure.children:
                if len(child.children)==0:
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
                                if len(sub_child.children)==0 and (child.name not in residues):
                                    warn('ERROR: All the residues are not specified')
                                    return None, None, None, None

        structure = new_compound_iter



    residues_applied_list = []
    Residue_orig_order_list = []
    for child in structure.children:
        if child.name not in  Residue_orig_order_list:
            Residue_orig_order_list.append(child.name)
    for res_reorder_iter in range(0,len(residues)):
        if residues[res_reorder_iter] not in  Residue_orig_order_list:
            text_to_print_1 = "All the residues were not used from the forcefield_names or forcefield_files "+ \
                              "string or dictionary.  There may be residues below other specified residues "+ \
                              "in the mbuild.Compound hierarchy.  If so, all the highest listed residues pass "+ \
                              "down the force fields through the hierarchy.  Alternatively, "+ \
                              "residues that are not in the structure may have been specified. "
            text_to_print_2 = "Note: This warning will appear if you are using the CHARMM pdb and psf writers "+ \
                              "2 boxes, and the boxes do not contain all the residues in each box."
            if boxes_for_simulation == 1:
                warn(text_to_print_1)
                return None, None, None, None
            if boxes_for_simulation == 2:
                warn(text_to_print_1 +text_to_print_2)

    if reorder_res_in_pdb_psf==False:
        residues = Residue_orig_order_list
    elif reorder_res_in_pdb_psf==True:
        print("Information: the output file are being reordered in via the residues list's sequence. ")
    else:
        warn("ERROR residues = Residue_orig_order_list or residues= residues not properly specified ")
        return None, None, None, None



    for i in range(0, len(residues)):
        children_in_iteration = False
        new_compound_iteration = mb.Compound()
        new_compound_iter.periodicity[0] = structure.periodicity[0]
        new_compound_iter.periodicity[1] = structure.periodicity[1]
        new_compound_iter.periodicity[2] = structure.periodicity[2]
        new_structure_iteration = pmd.Structure()
        new_structure_iteration.box = compound_box_infor.box
        for child in structure.children:
            if FF_data.get(child.name) is None:
                warn('ERROR: All residues are not specified in the force_field dictionary')
                return None, None, None, None

            if child.name == residues[i]:
                children_in_iteration = True
                new_compound_iteration.add(mb.compound.clone(child))

        if children_in_iteration == True:
            if Use_FF_files == True:
                FF_iteration = Forcefield(FF_data[residues[i]])
                residues_applied_list.append(residues[i])
            elif Use_FF_names == True:
                FF_iteration = Forcefield(name=FF_data[residues[i]])
                residues_applied_list.append(residues[i])


            new_structure_iteration = FF_iteration.apply(new_compound_iteration, residues=[residues[i]])
            new_structure = new_structure + new_structure_iteration

    if box != None:
        new_structure.box[0] = box_Ang[0]
        new_structure.box[1] = box_Ang[1]
        new_structure.box[2] = box_Ang[2]


    structure = new_structure

    return structure, coulomb14scaler_dict, LJ14scaler_dict, residues_applied_list