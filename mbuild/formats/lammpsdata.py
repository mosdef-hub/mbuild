from __future__ import division

__all__ = ['write_lammpsdata']


import numpy as np
from .hoomdxml import RB_to_OPLS

def write_lammpsdata(structure, filename, forcefield, box):
    """Output a LAMMPS data file.
    
    Note: Output supports 'real' units and 'full' atom style only.

    Parameters
    ----------
    structure : parmed.Structure
        Parmed structure object
    filename : str
        Path of the output file
    forcefield : str
        Name of the force field to be applied to the compound
    box : mb.Box
        Box information to save to data file
    """

    masses = [atom.mass for atom in structure.atoms]

    if forcefield:
        types_str = [atom.type for atom in structure.atoms]
        types_str_num = [int(atom.split('_')[1]) for atom in types_str]
        header = types_str[0].split('_')[0]
        all_types_num = list(set(types_str_num))
        all_types_num.sort()
        all_types_num = np.array([[num+1,atom_type] for num,atom_type in enumerate(all_types_num)])
        all_types = np.array([[num,header+'_'+str(atom_type)] for num,atom_type in all_types_num])
        type_dict_num = {pair[1]:pair[0] for pair in all_types_num}

        mass_dict = dict([(type_dict_num[atom_type],mass) for atom_type,mass in zip(types_str_num,masses)])
        types = [type_dict_num[atom_type] for atom_type in types_str_num]
    else:   
        types_str = [atom.name for atom in structure.atoms]
        all_types = list(set(types_str))
        all_types = np.array([[num+1,atom_type] for num,atom_type in enumerate(all_types)])
        type_dict = {pair[1]:int(pair[0]) for pair in all_types}
        mass_dict = dict([(type_dict[atom_type],mass) for atom_type,mass in zip(types_str,masses)])
        types = [type_dict[atom_type] for atom_type in types_str]

    xyz = np.array([[atom.xx,atom.xy,atom.xz] for atom in structure.atoms])
    charges = [atom.charge for atom in structure.atoms]

    bonds = [[bond.atom1.idx+1, bond.atom2.idx+1] for bond in structure.bonds]
    angles = [[angle.atom1.idx+1,
               angle.atom2.idx+1,
               angle.atom3.idx+1] for angle in structure.angles]
    dihedrals = [[dihedral.atom1.idx+1,
                  dihedral.atom2.idx+1,
                  dihedral.atom3.idx+1,
                  dihedral.atom4.idx+1] for dihedral in structure.rb_torsions]

    if bonds:
        if len(structure.bond_types) == 0:
            bond_types = np.ones(len(bonds),dtype=int)
        else:
            all_bond_types = dict(enumerate(set([(round(bond.type.k,3),
                                                  round(bond.type.req,3)) for bond in structure.bonds])))
            all_bond_types = {y:x+1 for x,y in all_bond_types.items()}
            bond_types = [all_bond_types[(round(bond.type.k,3),
                                          round(bond.type.req,3))] for bond in structure.bonds]

    if angles:
        all_angle_types = dict(enumerate(set([(round(angle.type.k,3),
                                               round(angle.type.theteq,3)) for angle in structure.angles])))
        all_angle_types = {y:x+1 for x,y in all_angle_types.items()}
        angle_types = [all_angle_types[(round(angle.type.k,3),
                                        round(angle.type.theteq,3))] for angle in structure.angles]

    if dihedrals:
        all_dihedral_types = dict(enumerate(set([(round(dihedral.type.c0,3),
                                                  round(dihedral.type.c1,3),
                                                  round(dihedral.type.c2,3),
                                                  round(dihedral.type.c3,3),
                                                  round(dihedral.type.c4,3),
                                                  round(dihedral.type.c5,3),
                                                  round(dihedral.type.scee,1),
                                                  round(dihedral.type.scnb,1)) for dihedral in structure.rb_torsions])))
        all_dihedral_types = {y:x+1 for x,y in all_dihedral_types.items()}
        dihedral_types = [all_dihedral_types[(round(dihedral.type.c0,3),
                                              round(dihedral.type.c1,3),
                                              round(dihedral.type.c2,3),
                                              round(dihedral.type.c3,3),
                                              round(dihedral.type.c4,3),
                                              round(dihedral.type.c5,3),
                                              round(dihedral.type.scee,1),
                                              round(dihedral.type.scnb,1))] for dihedral in structure.rb_torsions]

    with open(filename, 'w') as data:
        data.write(filename+' - created by mBuild\n\n')
        data.write('{:d} atoms\n'.format(len(structure.atoms)))
        data.write('{:d} bonds\n'.format(len(bonds)))
        data.write('{:d} angles\n'.format(len(angles)))
        data.write('{:d} dihedrals\n\n'.format(len(dihedrals)))

        data.write('{:d} atom types\n'.format(len(set(types))))
        if bonds:
            data.write('{:d} bond types\n'.format(len(set(bond_types))))
        if angles:
            data.write('{:d} angle types\n'.format(len(set(angle_types))))
        if dihedrals:
            data.write('{:d} dihedral types\n'.format(len(set(dihedral_types))))

        data.write('\n')
        # Box data
        for i,dim in enumerate(['x','y','z']):
            data.write('{0:.6f} {1:.6f} {2}lo {2}hi\n'.format(box.mins[i]*10.,box.maxs[i]*10.,dim))

        # Mass data
        type_dict_r = {int(pair[0]):pair[1] for pair in all_types}
        data.write('\nMasses\n\n')
        for atom_type,mass in mass_dict.items():
            data.write('{:d}\t{:.6f}\t# {}\n'.format(atom_type,mass,type_dict_r[atom_type]))

        if forcefield:
            # Pair coefficients
            epsilons = [atom.epsilon for atom in structure.atoms]
            sigmas = [atom.sigma for atom in structure.atoms]
            epsilon_dict = dict([(type_dict_num[atom_type],epsilon) for atom_type,epsilon in zip(types_str_num,epsilons)])
            sigma_dict = dict([(type_dict_num[atom_type],sigma) for atom_type,sigma in zip(types_str_num,sigmas)])
            data.write('\nPair Coeffs # lj\n\n')
            for idx,epsilon in epsilon_dict.items():
                data.write('{}\t{:.5f}\t{:.5f}\n'.format(idx,epsilon,sigma_dict[idx]))

            # Bond coefficients
            data.write('\nBond Coeffs # harmonic\n\n')
            for params,idx in all_bond_types.items():
                data.write('{}\t{}\t{}\n'.format(idx,*params))

            # Angle coefficients
            data.write('\nAngle Coeffs # harmonic\n\n')
            for params,idx in all_angle_types.items():
                data.write('{}\t{}\t{:.5f}\n'.format(idx,*params))

            # Dihedral coefficients
            data.write('\nDihedral Coeffs # opls\n\n')
            for params,idx in all_dihedral_types.items():
                opls_coeffs = RB_to_OPLS(params[0],
                                         params[1],
                                         params[2],
                                         params[3],
                                         params[4],
                                         params[5])
                data.write('{}\t{:.5f}\t{:.5f}\t{:.5f}\t{:.5f}\n'.format(idx,*opls_coeffs))

        # Atom data
        data.write('\nAtoms\n\n')
        for i,coords in enumerate(xyz):
            data.write('{:d}\t{:d}\t{:d}\t{:.6f}\t{:.6f}\t{:.6f}\t{:.6f}\n'.format(i+1,0,types[i],charges[i],*coords))

        # Bond data
        if bonds:
            data.write('\nBonds\n\n')
            for i,bond in enumerate(bonds):
                data.write('{:d}\t{:d}\t{:d}\t{:d}\n'.format(i+1,bond_types[i],bond[0],bond[1]))

        # Angle data
        if angles:
            data.write('\nAngles\n\n')
            for i,angle in enumerate(angles):
                data.write('{:d}\t{:d}\t{:d}\t{:d}\t{:d}\n'.format(i+1,angle_types[i],angle[0],angle[1],angle[2]))

        # Dihedral data
        if dihedrals:
            data.write('\nDihedrals\n\n')
            for i,dihedral in enumerate(dihedrals):
                data.write('{:d}\t{:d}\t{:d}\t{:d}\t{:d}\t{:d}\n'.format(i+1,dihedral_types[i],dihedral[0],dihedral[1],dihedral[2],dihedral[3]))
