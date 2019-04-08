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

    Some of this function has beed adopted from `mdtraj`'s support of the LAMMPSTRJ
    trajectory format. See https://github.com/mdtraj/mdtraj/blob/master/mdtraj/formats/lammpstrj.py for details.

    """

    if atom_style not in ['atomic', 'charge', 'molecular', 'full']:
        raise ValueError('Atom style "{}" is invalid or is not currently supported'.format(atom_style))

    xyz = np.array([[atom.xx,atom.xy,atom.xz] for atom in structure.atoms])

    forcefield = True
    if structure[0].type == '':
        forcefield = False

    # Internally use nm
    box = Box(lengths=np.array([0.1 * val for val in structure.box[0:3]]),
              angles=structure.box[3:6])

    if forcefield:
        types = [atom.type for atom in structure.atoms]
    else:
        types = [atom.name for atom in structure.atoms]

    unique_types = list(set(types))
    unique_types.sort(key=natural_sort)

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
            unique_bond_types = dict(enumerate(set([(round(bond.type.k,3),
                                                     round(bond.type.req,3)) for bond in structure.bonds])))
            unique_bond_types = OrderedDict([(y,x+1) for x,y in unique_bond_types.items()])
            bond_types = [unique_bond_types[(round(bond.type.k,3),
                                             round(bond.type.req,3))] for bond in structure.bonds]

    if angles:
        unique_angle_types = dict(enumerate(set([(round(angle.type.k,3),
                                                  round(angle.type.theteq,3)) for angle in structure.angles])))
        unique_angle_types = OrderedDict([(y,x+1) for x,y in unique_angle_types.items()])
        angle_types = [unique_angle_types[(round(angle.type.k,3),
                                           round(angle.type.theteq,3))] for angle in structure.angles]

    if dihedrals:
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

    with open(filename, 'w') as data:
        data.write(filename+' - created by mBuild\n\n')
        data.write('{:d} atoms\n'.format(len(structure.atoms)))
        if atom_style in ['full', 'molecular']:
            data.write('{:d} bonds\n'.format(len(bonds)))
            data.write('{:d} angles\n'.format(len(angles)))
            data.write('{:d} dihedrals\n\n'.format(len(dihedrals)))

        data.write('{:d} atom types\n'.format(len(set(types))))
        if atom_style in ['full', 'molecular']:
            if bonds:
                data.write('{:d} bond types\n'.format(len(set(bond_types))))
            if angles:
                data.write('{:d} angle types\n'.format(len(set(angle_types))))
            if dihedrals:
                data.write('{:d} dihedral types\n'.format(len(set(dihedral_types))))

        data.write('\n')
        # Box data
        if np.allclose(box.angles, np.array([90, 90, 90])):
            for i,dim in enumerate(['x','y','z']):
                data.write('{0:.6f} {1:.6f} {2}lo {2}hi\n'.format(
                    10.0 * box.mins[i],
                    10.0 * box.maxs[i],
                    dim))
        else:
            a, b, c = 10.0 * box.lengths
            alpha, beta, gamma = np.radians(box.angles)

            lx = a
            xy = b * np.cos(gamma)
            xz = c * np.cos(beta)
            ly = np.sqrt(b**2 - xy**2)
            yz = (b*c*np.cos(alpha) - xy*xz) / ly
            lz = np.sqrt(c**2 - xz**2 - yz**2)

            xlo, ylo, zlo = 10.0 * box.mins
            xhi = xlo + lx
            yhi = ylo + ly
            zhi = zlo + lz

            xlo_bound = xlo + np.min([0.0, xy, xz, xy+xz])
            xhi_bound = xhi + np.max([0.0, xy, xz, xy+xz])
            ylo_bound = ylo + np.min([0.0, yz])
            yhi_bound = yhi + np.max([0.0, yz])
            zlo_bound = zlo
            zhi_bound = zhi

            data.write('{0:.6f} {1:.6f} xlo xhi\n'.format(
                xlo_bound, xhi_bound))
            data.write('{0:.6f} {1:.6f} ylo yhi\n'.format(
                ylo_bound, yhi_bound))
            data.write('{0:.6f} {1:.6f} zlo zhi\n'.format(
                zlo_bound, zhi_bound))
            data.write('{0:.6f} {1:.6f} {2:6f} xy xz yz\n'.format(
                xy, xz, yz))

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
                data.write('\nDihedral Coeffs # opls\n\n')
                for params,idx in unique_dihedral_types.items():
                    opls_coeffs = RB_to_OPLS(params[0],
                                             params[1],
                                             params[2],
                                             params[3],
                                             params[4],
                                             params[5])
                    data.write('{}\t{:.5f}\t{:.5f}\t{:.5f}\t{:.5f}\n'.format(idx,*opls_coeffs))

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
