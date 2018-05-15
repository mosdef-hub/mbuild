from __future__ import division

from collections import OrderedDict

import numpy as np

from mbuild import Box
from mbuild.utils.conversion import RB_to_OPLS
from mbuild.utils.sorting import natural_sort

__all__ = ['write_lammpsdata']


def write_lammpsdata(structure, filename, atom_style='full'):
    """Output a LAMMPS data file.
    
    Outputs a LAMMPS data file in the 'full' atom style format. Assumes use
    of 'real' units. See http://lammps.sandia.gov/doc/atom_style.html for
    more information on atom styles.

    Parameters
    ----------
    structure : parmed.Structure
        ParmEd structure object
    filename : str
        Path of the output file
    atom_style: str
        Defines the style of atoms to be saved in a LAMMPS data file. The following atom
        styles are currently supported: 'full', 'atomic', 'charge', 'molecular'
        see http://lammps.sandia.gov/doc/atom_style.html for more
        information on atom styles.

    Notes
    -----
    See http://lammps.sandia.gov/doc/2001/data_format.html for a full description
    of the LAMMPS data format. Currently the following sections are supported (in
    addition to the header): *Masses*, *Nonbond Coeffs*, *Bond Coeffs*, *Angle
    Coeffs*, *Dihedral Coeffs*, *Atoms*, *Bonds*, *Angles*, *Dihedrals*

    """

    if atom_style not in ['atomic', 'charge', 'molecular', 'full']:
        raise ValueError('Atom style "{}" is invalid or is not currently supported'.format(atom_style))

    xyz = np.array([[atom.xx,atom.xy,atom.xz] for atom in structure.atoms])

    forcefield = True
    if structure[0].type == '':
        forcefield = False

    box = Box(lengths=np.array([structure.box[0], structure.box[1], structure.box[2]]))

    if forcefield:
        types = [atom.type for atom in structure.atoms]
    else:
        types = [atom.name for atom in structure.atoms]

    unique_types = list(set(types))
    unique_types.sort(key=natural_sort)

    charges = [atom.charge for atom in structure.atoms]

    
    # Lammps syntax depends on the functional form
    # Infer functional form based on the properties of the structure

    # Check angles
    if len(structure.urey_bradleys)> 0 :
        print("Urey bradley terms detected, will use angle_style charmm")
        use_urey_bradleys = True
    else:
        print("No urey bradley terms detected, will use angle_style harmonic")
        use_urey_bradleys = False

    # Check dihedrals
    if len(structure.rb_torsions) > 0:
        print("RB Torsions detected, will use dihedral_style opls")
        use_rb_torsions = True
    else:
        use_rb_torsions = False
    if len(structure.dihedrals) > 0:
        print("Charmm dihedrals detected, will use dihedral_style charmm")
        use_dihedrals = True
    else:
        use_dihedrals = False
    if use_rb_torsions and use_dihedrals:
        raise FoyerError("Multiple dihedral styles detected, check your "
                         "Forcefield XML and structure")

    # Check impropers
    for dihedral in structure.dihedrals:
        if dihedral.improper:
            raise FoyerError("Amber-style impropers are currently not supported")

    bonds = [[bond.atom1.idx+1, bond.atom2.idx+1] for bond in structure.bonds]
    angles = [[angle.atom1.idx+1,
               angle.atom2.idx+1,
               angle.atom3.idx+1] for angle in structure.angles]
    if use_rb_torsions:
        dihedrals = [[dihedral.atom1.idx+1,
                      dihedral.atom2.idx+1,
                      dihedral.atom3.idx+1,
                      dihedral.atom4.idx+1] for dihedral in structure.rb_torsions]
    elif use_dihedrals:
        dihedrals = [[dihedral.atom1.idx+1,
                      dihedral.atom2.idx+1,
                      dihedral.atom3.idx+1,
                      dihedral.atom4.idx+1] for dihedral in structure.dihedrals]
    else:
        dihedrals = []
    impropers = [[improper.atom1.idx+1,
                  improper.atom2.idx+1,
                  improper.atom3.idx+1,
                  improper.atom4.idx+1] for improper in structure.impropers]



    if bonds:
        if len(structure.bond_types) == 0:
            bond_types = np.ones(len(bonds),dtype=int)
        else:
            unique_bond_types = dict(enumerate(set([(round(bond.type.k,3),
                                                     round(bond.type.req,3)) for bond in structure.bonds])))
            unique_bond_types = OrderedDict([(y,x+1) for x,y in unique_bond_types.items()])
            bond_types = [unique_bond_types[(round(bond.type.k,3),
                                             round(bond.type.req,3))] for bond in structure.bonds]

    if angles:
        if use_urey_bradleys:
            charmm_angle_types = []
            for angle in structure.angles:
                for ub in structure.urey_bradleys:
                    if (angle.atom1, angle.atom3) == (ub.atom1, ub.atom2):
                        ub_k = ub.type.k
                        ub_req = ub.type.req
                    else:
                        ub_k = 0
                        ub_req = 0
                    charmm_angle_types.append((round(angle.type.k,3), 
                                               round(angle.type.theteq,3),
                                               round(ub_k, 3),
                                               round(ub_req, 3)))

            unique_angle_types = dict(enumerate(set(charmm_angle_types)))
            unique_angle_types = OrderedDict([(y,x+1) for x,y in unique_angle_types.items()])
            angle_types = [unique_angle_types[ub_info] for ub_info in charmm_angle_types]


        else:
            unique_angle_types = dict(enumerate(set([(round(angle.type.k,3),
                                                      round(angle.type.theteq,3)) for angle in structure.angles])))
            unique_angle_types = OrderedDict([(y,x+1) for x,y in unique_angle_types.items()])
            angle_types = [unique_angle_types[(round(angle.type.k,3),
                                               round(angle.type.theteq,3))] for angle in structure.angles]

    if dihedrals:
        if use_rb_torsions:
            unique_dihedral_types = dict(enumerate(set([(round(dihedral.type.c0,3),
                                                         round(dihedral.type.c1,3),
                                                         round(dihedral.type.c2,3),
                                                         round(dihedral.type.c3,3),
                                                         round(dihedral.type.c4,3),
                                                         round(dihedral.type.c5,3),
                                                         round(dihedral.type.scee,1),
                                                         round(dihedral.type.scnb,1)) for dihedral in structure.rb_torsions])))
            unique_dihedral_types = OrderedDict([(y,x+1) for x,y in unique_dihedral_types.items()])
            dihedral_types = [unique_dihedral_types[(round(dihedral.type.c0,3),
                                                     round(dihedral.type.c1,3),
                                                     round(dihedral.type.c2,3),
                                                     round(dihedral.type.c3,3),
                                                     round(dihedral.type.c4,3),
                                                     round(dihedral.type.c5,3),
                                                     round(dihedral.type.scee,1),
                                                     round(dihedral.type.scnb,1))] for dihedral in structure.rb_torsions]
        elif use_dihedrals:
            charmm_dihedrals = []
            for dihedral in structure.dihedrals:
                if not dihedral.improper:
                    weight = 1 / len(dihedral.type.list)
                    charmm_dihedrals.append((round(dihedral.type.phi_k,3),
                                             int(round(dihedral.type.per,0)),
                                             int(round(dihedral.type.phase,0)),
                                             round(weight, 1),
                                             round(dihedral.type.scee,1),
                                             round(dihedral.type.scnb,1)))

            unique_dihedral_types = dict(enumerate(set(charmm_dihedrals)))
            unique_dihedral_types = OrderedDict([(y,x+1) for x,y in unique_dihedral_types.items()])
            dihedral_types = [unique_dihedral_types[dihedral_info] for dihedral_info in charmm_dihedrals]
            
    if impropers:
            unique_improper_types = dict(enumerate(set([(round(improper.type.psi_k,3),
                                                         round(improper.type.psi_eq,3)) for improper in structure.impropers])))
            unique_improper_types = OrderedDict([(y,x+1) for x,y in unique_improper_types.items()])
            improper_types = [unique_improper_types[(round(improper.type.psi_k,3),
                                                     round(improper.type.psi_eq,3))] for improper in structure.impropers]

    with open(filename, 'w') as data:
        data.write(filename+' - created by mBuild\n\n')
        data.write('{:d} atoms\n'.format(len(structure.atoms)))
        if atom_style in ['full', 'molecular']:
            data.write('{:d} bonds\n'.format(len(bonds)))
            data.write('{:d} angles\n'.format(len(angles)))
            data.write('{:d} dihedrals\n'.format(len(dihedrals)))
            data.write('{:d} impropers\n\n'.format(len(impropers)))

        data.write('{:d} atom types\n'.format(len(set(types))))
        if atom_style in ['full', 'molecular']:
            if bonds:
                data.write('{:d} bond types\n'.format(len(set(bond_types))))
            if angles:
                data.write('{:d} angle types\n'.format(len(set(angle_types))))
            if dihedrals:
                data.write('{:d} dihedral types\n'.format(len(set(dihedral_types))))
            if impropers:
                data.write('{:d} improper types\n'.format(len(set(improper_types))))


        data.write('\n')
        # Box data
        for i,dim in enumerate(['x','y','z']):
            data.write('{0:.6f} {1:.6f} {2}lo {2}hi\n'.format(box.mins[i],box.maxs[i],dim))

        # Mass data
        masses = [atom.mass for atom in structure.atoms]
        mass_dict = dict([(unique_types.index(atom_type)+1,mass) for atom_type,mass in zip(types,masses)])

        data.write('\nMasses\n\n')
        for atom_type,mass in mass_dict.items():
            data.write('{:d}\t{:.6f}\t# {}\n'.format(atom_type,mass,unique_types[atom_type-1]))

        if forcefield:
            # Pair coefficients
            epsilons = [atom.epsilon for atom in structure.atoms]
            sigmas = [atom.sigma for atom in structure.atoms]
            epsilon_dict = dict([(unique_types.index(atom_type)+1,epsilon) for atom_type,epsilon in zip(types,epsilons)])
            sigma_dict = dict([(unique_types.index(atom_type)+1,sigma) for atom_type,sigma in zip(types,sigmas)])
            data.write('\nPair Coeffs # lj\n\n')
            for idx,epsilon in epsilon_dict.items():
                data.write('{}\t{:.5f}\t{:.5f}\n'.format(idx,epsilon,sigma_dict[idx]))

            # Bond coefficients
            if bonds:
                data.write('\nBond Coeffs # harmonic\n\n')
                for params,idx in unique_bond_types.items():
                    data.write('{}\t{}\t{}\n'.format(idx,*params))

            # Angle coefficients
            if angles:
                data.write('\nAngle Coeffs # harmonic\n\n')
                for params,idx in unique_angle_types.items():
                    data.write('{}\t{}\t{:.5f}\n'.format(idx,*params))

            # Dihedral coefficients
            if dihedrals:
                if use_rb_torsions:
                    data.write('\nDihedral Coeffs # opls\n\n')
                    for params,idx in unique_dihedral_types.items():
                        opls_coeffs = RB_to_OPLS(params[0],
                                                 params[1],
                                                 params[2],
                                                 params[3],
                                                 params[4],
                                                 params[5])
                        data.write('{}\t{:.5f}\t{:.5f}\t{:.5f}\t{:.5f}\n'.format(idx,*opls_coeffs))
                elif use_dihedrals:
                    data.write('\nDihedral Coeffs # charmm\n\n')
                    for params, idx in unique_dihedral_types.items():
                        data.write('{}\t{:.5f}\t{:d}\t{:d}\t{:.5f}\n'.format(idx, *params))

            # Improper coefficients
            if impropers:
                data.write('\nImproper Coeffs # harmonic\n\n')
                for params,idx in unique_improper_types.items():
                    data.write('{}\t{:.5f}\t{:.5f}\n'.format(idx, *params))

        # Atom data
        data.write('\nAtoms\n\n')
        if atom_style == 'atomic':
            atom_line = '{index:d}\t{type_index:d}\t{x:.6f}\t{y:.6f}\t{z:.6f}\n'
        elif atom_style == 'charge':
            atom_line = '{index:d}\t{type_index:d}\t{charge:.6f}\t{x:.6f}\t{y:.6f}\t{z:.6f}\n'
        elif atom_style == 'molecular':
            atom_line = '{index:d}\t{zero:d}\t{type_index:d}\t{x:.6f}\t{y:.6f}\t{z:.6f}\n'
        elif atom_style == 'full':
            atom_line ='{index:d}\t{zero:d}\t{type_index:d}\t{charge:.6f}\t{x:.6f}\t{y:.6f}\t{z:.6f}\n'

        for i,coords in enumerate(xyz):
            data.write(atom_line.format(
                index=i+1,type_index=unique_types.index(types[i])+1,
                zero=0,charge=charges[i],
                x=coords[0],y=coords[1],z=coords[2]))

        if atom_style in ['full', 'molecular']:
            # Bond data
            if bonds:
                data.write('\nBonds\n\n')
                for i,bond in enumerate(bonds):
                    data.write('{:d}\t{:d}\t{:d}\t{:d}\n'.format(
                        i+1,bond_types[i],bond[0],bond[1]))

            # Angle data
            if angles:
                data.write('\nAngles\n\n')
                for i,angle in enumerate(angles):
                    data.write('{:d}\t{:d}\t{:d}\t{:d}\t{:d}\n'.format(
                        i+1,angle_types[i],angle[0],angle[1],angle[2]))

            # Dihedral data
            if dihedrals:
                data.write('\nDihedrals\n\n')
                for i,dihedral in enumerate(dihedrals):
                    data.write('{:d}\t{:d}\t{:d}\t{:d}\t{:d}\t{:d}\n'.format(
                        i+1,dihedral_types[i],dihedral[0],
                        dihedral[1],dihedral[2],dihedral[3]))
            # Dihedral data
            if impropers:
                data.write('\nImpropers\n\n')
                for i,improper in enumerate(impropers):
                    data.write('{:d}\t{:d}\t{:d}\t{:d}\t{:d}\t{:d}\n'.format(
                        i+1,improper_types[i],improper[0],
                        improper[1],improper[2],improper[3]))
