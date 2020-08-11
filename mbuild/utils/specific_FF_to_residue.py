from foyer import Forcefield
import parmed as pmd
import mbuild as mb
import xml.dom.minidom
from warnings import warn

def Specific_FF_to_residue(structure , forcefield_files= None, residues= None, reorder_res_in_pdb_psf= False):

    #input:
        # structure =  compound structure
        #  forcefield_files =  dictionary of residues to force fields.
            # Ex: {'Water' : 'oplsaa.xml', 'OCT': 'trappe-ua.xml}
        # residues = list of residues
        # note FFs and residues must be in sequenctial order

    #Returns:
        # structure = parmed with applied force field



    if  forcefield_files== None:
        print('please enter the forcefield_files in the Specific_FF_to_residue function')
    if forcefield_files != None and not isinstance(forcefield_files, dict):
        return warn('The force field file (forcefield_files) is not a dictionary. Please enter a dictionary'+
                      "with all the residues specified to a force field"+
                     "-> Ex: {'Water' : 'oplsaa.xml', 'OCT': 'trappe-ua.xml'}, "+
                     "Note: the file path must be specified the force field file")

    if residues == None:
        print('please enter the residues in the Specific_FF_to_residue function')
    if reorder_res_in_pdb_psf == None:
        print('please enter the reorder_res_in_pdb_psf in the Specific_FF_to_residue function')


    # select 1,4-interatactions for each residue
    #forcefield_files_keys = list(forcefield_files.keys())
    #forcefield_files_values = forcefield_files.values()

    forcefield_files_keys_list = []

    for res in forcefield_files.keys():
        forcefield_files_keys_list.append(res)

    coulomb14scaler_dict = {}
    LJ14scaler_dict = {}
    for j in range(0, len(forcefield_files_keys_list)):
        residue_iteration = forcefield_files_keys_list[j]
        FF_for_residue_iteration = forcefield_files[residue_iteration]
        read_xlm_iteration = xml.dom.minidom.parse(FF_for_residue_iteration)
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

            #new_compound_iter.name = str(residues[i])
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
                                    #new_compound_iter.add(mb.compound.clone(child))


            structure = new_compound_iter

        """
        for child in structure.children:
            if child.name not in residues:
                print('All the residues are not specified')
                return warn('All the residues are not specified')
        """

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
            #print("filled_pore.child = " + str(child))
            if forcefield_files.get(child.name) == None:
                print('All residues are not specified in the force_field dictionary')
                return warn('All residues are not specified in the force_field dictionary')

            if child.name == residues[i]:
                children_in_iteration = True
                new_compound_iteration.add(mb.compound.clone(child))

        if children_in_iteration == True:
            FF_iteration = Forcefield(forcefield_files[residues[i]])

            new_structure_iteration = FF_iteration.apply(new_compound_iteration, residues=[residues[i]])
            new_structure = new_structure + new_structure_iteration
            #print('new_structure = ' + str(new_structure))

    structure = new_structure
    return structure, coulomb14scaler_dict, LJ14scaler_dict

