from mbuild.prototype import Prototype
from mbuild.treeview import TreeView

__author__ = 'sallai'
from mbuild.compound import *
import os.path
import re


class Lammps(Compound):

    def __init__(self, path, ctx={}, cwd="", verbose=False):
        """Reads a LAMMPS data file

            *** Only works for directives delimited by blank lines ***

            Currently supports the following directives:
                Masses
                Pair Coeffs (must be mix geometry)
                Bond Coeffs (must be harmonic)
                Angle Coeffs (must be harmonic)
                Dihedral Coeffs (must be OPLS)
                Atoms
                Bonds
                Angles
                Dihedrals

            TODO:
                -handling for comments
                -handling for directives not delimited by blank lines
                -allow specification of forcefield styles

            Args:
                data_file (str): name of LAMMPS data file to read in
            Returns:
                lmp_data (dict):
                    'xyz': xyz (numpy.ndarray)
                    'types': types (numpy.ndarray)
                    'masses': masses (numpy.ndarray)
                    'charges': charges (numpy.ndarray)
                    'bonds': bonds (numpy.ndarray)
                    'angles': angles (numpy.ndarray)
                    'dihedrals': dihedrals (numpy.ndarray)
                    'pair_types': pair_types (dict)
                    'bond_types': bond_types (dict)
                    'angle_types': angle_types (dict)
                    'dihedral_types': dihedral_type (dict)

                box (numpy.ndarray): box dimensions
            """
        super(Lammps, self).__init__(ctx=ctx)

        data_file = os.path.join(cwd, path)

        bonds = np.empty(shape=(0, 3), dtype='int')
        angles = np.empty(shape=(0, 4), dtype='int')
        dihedrals = np.empty(shape=(0, 5), dtype='int')
        impropers = np.empty(shape=(0, 5), dtype='int')

        pair_types = dict()
        bond_types = dict()
        angle_types = dict()
        dihedral_types = dict()
        improper_types = dict()

        print "Reading '" + data_file + "'"
        with open(data_file, 'r') as f:
            data_lines = f.readlines()

        # TODO: improve robustness of xlo regex
        directives = re.compile(r"""
            ((?P<n_atoms>\s*\d+\s+atoms)
            |
            (?P<n_bonds>\s*\d+\s+bonds)
            |
            (?P<n_angles>\s*\d+\s+angles)
            |
            (?P<n_dihedrals>\s*\d+\s+dihedrals)
            |
            (?P<box>.+xlo)
            |
            (?P<Masses>\s*Masses)
            |
            (?P<PairCoeffs>\s*Pair\sCoeffs)
            |
            (?P<BondCoeffs>\s*Bond\sCoeffs)
            |
            (?P<AngleCoeffs>\s*Angle\sCoeffs)
            |
            (?P<DihedralCoeffs>\s*Dihedral\sCoeffs)
            |
            (?P<ImproperCoeffs>\s*Improper\sCoeffs)
            |
            (?P<Atoms>\s*Atoms)
            |
            (?P<Bonds>\s*Bonds)
            |
            (?P<Angles>\s*Angles)
            |
            (?P<Dihedrals>\s*Dihedrals)
            |
            (?P<Impropers>\s*Impropers))
            """, re.VERBOSE)

        i = 0
        while i < len(data_lines):
            match = directives.match(data_lines[i])
            if match:
                if verbose:
                    print(match.groups())

                elif match.group('n_atoms'):
                    fields = data_lines.pop(i).split()
                    n_atoms = int(fields[0])
                    xyz = np.empty(shape=(n_atoms, 3))
                    types = np.empty(shape=(n_atoms), dtype='int')
                    masses = np.empty(shape=(n_atoms))
                    charges = np.empty(shape=(n_atoms))

                elif match.group('n_bonds'):
                    fields = data_lines.pop(i).split()
                    bonds = np.empty(shape=(float(fields[0]), 3), dtype='int')

                elif match.group('n_angles'):
                    fields = data_lines.pop(i).split()
                    angles = np.empty(shape=(float(fields[0]), 4), dtype='int')

                elif match.group('n_dihedrals'):
                    fields = data_lines.pop(i).split()
                    dihedrals = np.empty(shape=(float(fields[0]), 5), dtype='int')

                elif match.group('box'):
                    dims = np.zeros(shape=(3, 2))
                    for j in range(3):
                        fields = map(float, data_lines.pop(i).split()[:2])
                        dims[j, 0] = fields[0]
                        dims[j, 1] = fields[1]
                    box_lo = dims[:, 0]
                    box_hi = dims[:, 1]

                elif match.group('Masses'):
                    if verbose:
                        print 'Parsing Masses...'
                    data_lines.pop(i)
                    data_lines.pop(i)

                    mass_dict = dict()  # type:mass
                    #     not end of file         not blank line
                    while i < len(data_lines) and data_lines[i].strip():
                        fields = data_lines.pop(i).split()
                        mass_dict[int(fields[0])] = float(fields[1])

                elif match.group('Atoms'):
                    if verbose:
                        print 'Parsing Atoms...'
                    data_lines.pop(i)
                    data_lines.pop(i)

                    while i < len(data_lines) and data_lines[i].strip():
                        fields = data_lines.pop(i).split()
                        if len(fields) == 7:
                            a_id = int(fields[0])
                            types[a_id - 1] = int(fields[2])
                            masses[a_id - 1] = mass_dict[int(fields[2])]
                            charges[a_id - 1] = float(fields[3])
                            xyz[a_id - 1] = np.array([float(fields[4]),
                                                 float(fields[5]),
                                                 float(fields[6])])

                        # non-official file format
                        if len(fields) == 8:
                            a_id = int(fields[0])
                            types[a_id - 1] = int(fields[1])
                            masses[a_id - 1] = mass_dict[int(fields[1])]
                            charges[a_id - 1] = 0.0
                            xyz[a_id - 1] = np.array([float(fields[2]),
                                                 float(fields[3]),
                                                 float(fields[4])])


                elif match.group('PairCoeffs'):
                    if verbose:
                        print 'Parsing Pair Coeffs...'
                    data_lines.pop(i)
                    data_lines.pop(i)

                    while i < len(data_lines) and data_lines[i].strip():
                        fields = data_lines.pop(i).split()
                        pair_types[int(fields[0])] = (float(fields[1]),
                                                      float(fields[2]))
                elif match.group('BondCoeffs'):
                    if verbose:
                        print 'Parsing Bond Coeffs...'
                    data_lines.pop(i)
                    data_lines.pop(i)

                    while i < len(data_lines) and data_lines[i].strip():
                        fields = map(float, data_lines.pop(i).split())
                        bond_types[int(fields[0])] = fields[1:]

                elif match.group('AngleCoeffs'):
                    if verbose:
                        print 'Parsing Angle Coeffs...'
                    data_lines.pop(i)
                    data_lines.pop(i)

                    while i < len(data_lines) and data_lines[i].strip():
                        fields = map(float, data_lines.pop(i).split())
                        angle_types[int(fields[0])] = fields[1:]

                elif match.group('DihedralCoeffs'):
                    if verbose:
                        print 'Parsing Dihedral Coeffs...'
                    data_lines.pop(i)
                    data_lines.pop(i)

                    while i < len(data_lines) and data_lines[i].strip():
                        fields = map(float, data_lines.pop(i).split())
                        dihedral_types[int(fields[0])] = fields[1:]

                elif match.group('ImproperCoeffs'):
                    if verbose:
                        print 'Parsing Improper Coeffs...'
                    data_lines.pop(i)
                    data_lines.pop(i)

                    while i < len(data_lines) and data_lines[i].strip():
                        fields = map(float, data_lines.pop(i).split())
                        improper_types[int(fields[0])] = fields[1:]

                elif match.group('Bonds'):
                    if verbose:
                        print 'Parsing Bonds...'
                    data_lines.pop(i)
                    data_lines.pop(i)

                    while i < len(data_lines) and data_lines[i].strip():
                        fields = map(int, data_lines.pop(i).split())
                        bonds[fields[0] - 1] = fields[1:]

                elif match.group('Angles'):
                    if verbose:
                        print 'Parsing Angles...'
                    data_lines.pop(i)
                    data_lines.pop(i)

                    while i < len(data_lines) and data_lines[i].strip():
                        fields = map(int, data_lines.pop(i).split())
                        angles[fields[0] - 1] = fields[1:]

                elif match.group('Dihedrals'):
                    if verbose:
                        print 'Parsing Dihedrals...'
                    data_lines.pop(i)
                    data_lines.pop(i)

                    while i < len(data_lines) and data_lines[i].strip():
                        fields = map(int, data_lines.pop(i).split())
                        dihedrals[fields[0] - 1] = fields[1:]

                elif match.group('Impropers'):
                    if verbose:
                        print 'Parsing Impropers...'
                    data_lines.pop(i)
                    data_lines.pop(i)

                    while i < len(data_lines) and data_lines[i].strip():
                        fields = map(int, data_lines.pop(i).split())
                        dihedrals[fields[0] - 1] = fields[1:]

                else:
                    i += 1
            else:
                i += 1

        lmp_data = {'xyz': xyz,
                    'types': types,
                    'masses': masses,
                    'charges': charges,
                    'bonds': bonds,
                    'angles': angles,
                    'dihedrals': dihedrals,
                    'impropers': impropers,
                    'pair_types': pair_types,
                    'bond_types': bond_types,
                    'angle_types': angle_types,
                    'dihedral_types': dihedral_types,
                    'improper_types': improper_types
                    }

        # copy over everything relevant to the compound's data structures

        self.bounds = box_hi - box_lo
        xyz -= box_lo

        atoms_by_index = dict()
        for atom_idx, coords in enumerate(xyz):
            kind=str(types[atom_idx])
            a = Atom(kind=kind, pos=tuple(coords))
            a.charge = charges[atom_idx]
            atoms_by_index[atom_idx] = a
            self.add(a)
            Prototype(kind, mass=masses[atom_idx], pair_coeffs=pair_types[types[atom_idx]])

        # bonds_by_index = dict()
        for bond_idx, params in enumerate(bonds):
            lammps_bond_type = params[0]
            atom1_idx = params[1]-1
            atom2_idx = params[2]-1
            kind = "{0}-{1}".format(atom1_idx,atom2_idx)
            b = Bond(atoms_by_index[atom1_idx], atoms_by_index[atom2_idx], kind=kind)
            # bonds_by_index[bond_idx] = b
            self.add(b)
            bond_coeffs=bond_types[lammps_bond_type]
            Prototype(kind, bond_coeffs=bond_coeffs, lammps_bond_type=lammps_bond_type)

        # angles_by_index = dict()
        for angle_idx, params in enumerate(angles):
            kind = "{0}-{1}-{2}".format(types[params[1]-1],types[params[2]-1],types[params[3]-1])
            a = Angle(atoms_by_index[params[1]-1], atoms_by_index[params[2]-1], atoms_by_index[params[3]-1], kind=kind)
            # angles_by_index[angle_idx] = a
            self.add(a)
            Prototype(kind, angle_coeffs=angle_types[params[0]], lammps_angle_type=params[0])

        # dihedrals_by_index = dict()
        for dihedral_idx, params in enumerate(dihedrals):
            lammps_dihedral_type = params[0]
            atom1_idx = params[1]-1
            atom2_idx = params[2]-1
            atom3_idx = params[3]-1
            atom4_idx = params[4]-1

            print params

            kind = "{0}-{1}-{2}-{3}".format(types[params[1]-1],types[params[2]-1],types[params[3]-1],types[params[4]-1])
            d = Dihedral(atoms_by_index[params[1]-1], atoms_by_index[params[2]-1], atoms_by_index[params[3]-1], atoms_by_index[params[4]-1], kind=kind)
            # hidedrals_by_index[dihedral_idx] = d
            self.add(d)
            Prototype(kind, dihedral_coeffs=dihedral_types[params[0]], lammps_dihedral_type=params[0])

        # TODO: impropers
        # # impropers_by_index = dict()
        # for dihedral_idx, params in enumerate(dihedrals):
        #     kind = "{0}-{1}-{2}-{3}".format(types[params[1]-1],types[params[2]-1],types[params[3]-1],types[params[4]-1])
        #     d = Dihedral(atoms_by_index[params[1]-1], atoms_by_index[params[2]-1], atoms_by_index[params[3]-1], atoms_by_index[params[4]-1], kind=kind)
        #     # hidedrals_by_index[dihedral_idx] = d
        #     self.add(d)



    @staticmethod
    def save(compound, path, cwd="", sys_name="system", print_ports=False):
        """Write gbb to LAMMPS data file

        Args:
            compound (Compound): compound to write
            path (str): name of output file
            cwd (str): working directory
            sys_name (str): name printed at top of data file
        """

        fn = os.path.join(cwd, path)

        with open(fn, 'w') as f:
            f.write(sys_name + '\n')
            f.write('\n')

            n_atoms = len(compound.atoms)
            n_bonds = len(compound.bonds)
            n_angles = len(compound.angles)
            n_dihedrals = len(compound.dihedrals)

            f.write(str(n_atoms) + ' atoms\n')
            f.write(str(n_bonds) + ' bonds\n')
            f.write(str(n_angles) + ' angles\n')
            f.write(str(n_dihedrals) + ' dihedrals\n')
            f.write('\n')

            a_types, bond_types, ang_types, dih_types = compound.unique_types()

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

            f.write('0.0 %8.4f xlo xhi\n' % (compound.bounds[0]))
            f.write('0.0 %8.4f ylo yhi\n' % (compound.bounds[1]))
            f.write('0.0 %8.4f zlo zhi\n' % (compound.bounds[2]))

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

            for i, atom in enumerate(compound.atoms):
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
                for i, bond in enumerate(compound.bonds):
                    n_i = bond.atom1.id_num
                    n_j = bond.atom2.id_num
                    f.write('%d %d %d\n' % (i+1, n_i, n_j))

            if n_angles > 0:
                f.write('\n')
                f.write('Angles\n')
                f.write('\n')
                for i, angle in enumerate(compound.angles):
                    n_i = angle.atom1.id_num
                    n_j = angle.atom2.id_num
                    n_k = angle.atom3.id_num
                    f.write('%d %d %d %d\n' % (i+1, n_i, n_j, n_k))

            if n_dihedrals > 0:
                f.write('\n')
                f.write('Dihedrals\n')
                f.write('\n')
                for i, dihedral in enumerate(compound.dihedrals):
                    n_i = dihedral.atom1.id_num
                    n_j = dihedral.atom2.id_num
                    n_k = dihedral.atom3.id_num
                    n_l = dihedral.atom4.id_num
                    f.write('%d %d %d %d %d\n' % (i+1, n_i, n_j, n_k, n_l))


if __name__ == "__main__":
    m = Lammps("bmim.lmp", cwd="..\\lammps_testing")
    print [a for a in m.atoms()]

    print m.boundingbox()

    TreeView(m, bonds=False, angles=True).show()
