import numpy as np
import warnings
import re
import pdb


def write_lammps_data(mm, file_name='data.system', sys_name='system'):
    """Write gbb to LAMMPS data file

    Args:
        mm (MoleculeModel): molecule model to write
        box (numpy.ndarray): box dimensions
        file_name (str): name of output file
        sys_name (str): name printed at top of data file
    """
    with open(file_name, 'w') as f:
        f.write(sys_name + '\n')
        f.write('\n')

        n_atoms = len(mm.atoms)
        n_bonds = len(mm.bonds)
        n_angles = len(mm.angles)
        n_dihedrals = len(mm.dihedrals)

        f.write(str(n_atoms) + ' atoms\n')
        f.write(str(n_bonds) + ' bonds\n')
        f.write(str(n_angles) + ' angles\n')
        f.write(str(n_dihedrals) + ' dihedrals\n')
        f.write('\n')

        a_types, bond_types, ang_types, dih_types = mm.unique_types()

        # assign number to each unique atomtype
        numbered_a_types = dict(zip(a_types.keys(), range(1, len(a_types.keys()) + 1)))

        f.write(str(len(a_types)) + ' atom types\n')
        if n_bonds > 0:
            f.write(str(len(bond_types)) + ' bond types\n')
        if n_angles > 0:
            f.write(str(len(ang_types)) + ' angle types\n')
        if n_dihedrals > 0:
            f.write(str(len(dih_types)) + ' dihedral types\n')
        f.write('\n')

        f.write('0.0 %8.4f xlo xhi\n' % (mm.bounds[0]))
        f.write('0.0 %8.4f ylo yhi\n' % (mm.bounds[1]))
        f.write('0.0 %8.4f zlo zhi\n' % (mm.bounds[2]))

        f.write('\n')
        f.write('Masses\n')
        f.write('\n')

        # find unique masses and corresponding atomtypes
        masses = set()
        bogus_masses = [1.0] * len(a_types)  # need to properly implement mass
        for a_type, mass in zip(numbered_a_types.values(), bogus_masses):
            # really there should a lookup for mass based on type here
            masses.add((a_type, mass))
        for mass in sorted(masses):
            f.write(" ".join(map(str, mass)) + '\n')

        f.write('\n')
        f.write('Atoms\n')
        f.write('\n')

        for i, atom in enumerate(mm.atoms):
            atom.id_num = i
            f.write('%-6d %-6d %-6d %5.3f %8.3f %8.3f %8.3f\n'
                % (i+1,
                   1,  # TODO: molecule numbering
                   numbered_a_types[atom.kind],
                   0.0, # TODO: charge lookups
                   atom.pos[0],
                   atom.pos[1],
                   atom.pos[2]))

        if n_bonds > 0:
            f.write('\n')
            f.write('Bonds\n')
            f.write('\n')
            for i, bond in enumerate(mm.bonds):
                n_i = bond.atom1.id_num
                n_j = bond.atom2.id_num
                f.write('%d %d %d\n' % (i+1, n_i, n_j))

        if n_angles > 0:
            f.write('\n')
            f.write('Angles\n')
            f.write('\n')
            for i, angle in enumerate(mm.angles):
                n_i = angle.atom1.id_num
                n_j = angle.atom2.id_num
                n_k = angle.atom3.id_num
                f.write('%d %d %d %d\n' % (i+1, n_i, n_j, n_k))

        if n_dihedrals > 0:
            f.write('\n')
            f.write('Dihedrals\n')
            f.write('\n')
            for i, dihedral in enumerate(mm.dihedrals):
                n_i = dihedral.atom1.id_num
                n_j = dihedral.atom2.id_num
                n_k = dihedral.atom3.id_num
                n_l = dihedral.atom4.id_num
                f.write('%d %d %d %d %d\n' % (i+1, n_i, n_j, n_k, n_l))

    print "Wrote file '" + file_name + "'"

