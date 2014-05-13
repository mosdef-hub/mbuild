import pdb
import time
import os.path
import re

import mbuild.unit as units
from mbuild.prototype import Prototype
from mbuild.treeview import TreeView
from mbuild.compound import *


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

        self.periodicity = box_hi - box_lo
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
    def save(compound, ff, path, cwd="", unit_set='real', sys_name="system", print_ports=False):
        """Write gbb to LAMMPS data file

        Args:
            compound (Compound): compound to write
            path (str): name of output file
            cwd (str): working directory
            sys_name (str): name printed at top of data file
        """
        data_file = os.path.join(cwd, path)

        RAD = units.radians
        DEGREE = units.degrees
        if unit_set == 'real':
            DIST = units.angstroms
            VEL = units.angstroms / units.femtosecond
            ENERGY = units.kilocalorie / units.mole
            MASS = units.grams / units.mole
            CHARGE = units.elementary_charge
            MOLE = units.mole
        else:
            raise Exception("Unsupported unit set specified: {0}".format(unit_set))

        # Containers for lines which are ultimately written to output files.
        mass_list = list()
        mass_list.append('\n')
        mass_list.append('Masses\n')
        mass_list.append('\n')

        pair_coeff_list = list()
        pair_coeff_list.append('\n')
        pair_coeff_list.append('Pair Coeffs\n')
        pair_coeff_list.append('\n')

        bond_coeffs = list()
        bond_coeffs.append('\n')
        bond_coeffs.append('Bond Coeffs\n')
        bond_coeffs.append('\n')

        angle_coeffs = list()
        angle_coeffs.append('\n')
        angle_coeffs.append('Angle Coeffs\n')
        angle_coeffs.append('\n')

        dihedral_coeffs = list()
        dihedral_coeffs.append('\n')
        dihedral_coeffs.append('Dihedral Coeffs\n')
        dihedral_coeffs.append('\n')

        improper_coeffs = list()
        improper_coeffs.append('\n')
        improper_coeffs.append('Improper Coeffs\n')
        improper_coeffs.append('\n')

        atom_list = list()
        atom_list.append('\n')
        atom_list.append('Atoms\n')
        atom_list.append('\n')

        bond_list = list()
        bond_list.append('\n')
        bond_list.append('Bonds\n')
        bond_list.append('\n')

        angle_list = list()
        angle_list.append('\n')
        angle_list.append('Angles\n')
        angle_list.append('\n')

        dihedral_list = list()
        dihedral_list.append('\n')
        dihedral_list.append('Dihedrals\n')
        dihedral_list.append('\n')

        improper_list = list()
        improper_list.append('\n')
        improper_list.append('Impropers\n')
        improper_list.append('\n')

        # dicts for type information
        atom_type_dict = dict()
        a_type_i = 1  # counter for atom types

        bond_style = set()
        bond_type_dict = dict()
        b_type_i = 1  # counter for bond types

        angle_style = set()
        angle_type_dict = dict()
        ang_type_i = 1

        dihedral_style = set()
        dihedral_type_dict = dict()
        dih_type_i = 1

        improper_style = set()
        improper_type_dict = dict()
        imp_type_i = 1

        for i, atom in enumerate(compound.atoms()):
            x_min = np.inf
            y_min = np.inf
            z_min = np.inf
            x_max = -np.inf
            y_max = -np.inf
            z_max = -np.inf
            # type, mass and pair coeffs
            if atom.kind != 'G':
                if atom.kind not in atom_type_dict:
                    atom_type_dict[atom.kind] = a_type_i
                    a_type_i += 1

                    atom_type = ff.atom_types[atom.kind]
                    mass = atom_type.mass.in_units_of(MASS)._value
                    epsilon = atom_type.epsilon.in_units_of(ENERGY)._value
                    sigma = atom_type.sigma.in_units_of(DIST)._value
                    # mass
                    mass_list.append('{0:d} {1:8.4f}\n'.format(a_type_i, mass))
                    # pair coeff
                    pair_coeff_list.append('{0:d} {1:8.4f} {2:8.4f}\n'.format(
                            a_type_i, epsilon, sigma))

                x_coord = atom.pos[0]
                y_coord = atom.pos[1]
                z_coord = atom.pos[2]
                # box dimensions
                if x_coord < x_min:
                    x_min = x_coord
                if x_coord > x_max:
                    x_max = x_coord
                if y_coord < y_min:
                    y_min = y_coord
                if y_coord > y_max:
                    y_max = y_coord
                if z_coord < z_min:
                    z_min = z_coord
                if z_coord > z_max:
                    z_max = z_coord

                # atom
                atom.id_num = i
                atom_list.append('{0:-6d} {1:-6d} {2:-6d} {3:5.8f} {4:8.5f} {5:8.5f} {6:8.5f}\n'.format(
                        i + 1, 1, atom_type_dict[atom.kind],
                        atom.charge,
                        x_coord, y_coord, z_coord))

        for i, bond in enumerate(compound.bonds):
            if bond.kind not in bond_type_dict:
                bond_type_dict[bond.kind] = b_type_i
                pair = tuple(bond.kind.split('-'))

                r, k = ff.bond_types[pair]
                # NOTE: k includes the factor of 0.5 for harmonic in LAMMPS
                bond_coeffs.append('{0:d} {1:18.8f} {2:18.8f}\n'.format(
                        b_type_i,
                        0.5 * k.in_units_of(ENERGY / (DIST * DIST))._value,
                        r.in_units_of(DIST)._value))
                b_type_i += 1
            bond_list.append('{0:-6d} {1:6d} {2:6d} {3:6d}\n'.format(
                    i + 1,
                    bond_type_dict[bond.kind],
                    bond.atom1.id_num,
                    bond.atom2.id_num))
        """
                # angle types
                print "Writing angles..."
                for j, angle in enumerate(mol_type.angleForceSet.itervalues()):
                    atom1 = mol_type.moleculeSet[0]._atoms[angle.atom1 - 1]
                    atomtype1 = atom1.bondtype
                    atom2 = mol_type.moleculeSet[0]._atoms[angle.atom2 - 1]
                    atomtype2 = atom2.bondtype
                    atom3 = mol_type.moleculeSet[0]._atoms[angle.atom3 - 1]
                    atomtype3 = atom3.bondtype

                    if isinstance(angle, Angle):
                        style = 'harmonic'
                        temp = AngleType(atomtype1, atomtype2, atomtype3,
                                angle.theta, angle.k)
                        # NOTE: k includes the factor of 0.5 for harmonic in LAMMPS
                        if temp not in angle_type_dict:
                            angle_type_dict[temp] = ang_type_i
                            angle_coeffs.append('{0:d} {1} {2:18.8f} {3:18.8f}\n'.format(
                                    ang_type_i,
                                    style, 
                                    0.5 * angle.k.in_units_of(self.ENERGY / (self.RAD*self.RAD))._value,
                                    angle.theta.in_units_of(self.DEGREE)._value))
                            ang_type_i += 1
                    else:
                        warn("Found unimplemented angle type for LAMMPS!")
                        continue

                    angle_list.append('{0:-6d} {1:6d} {2:6d} {3:6d} {4:6d}\n'.format(
                            i + j + 1,
                            angle_type_dict[temp],
                            angle.atom1 + offset,
                            angle.atom2 + offset,
                            angle.atom3 + offset))

                    angle_style.add(style)
                    if len(angle_style) > 1:
                        warn("More than one angle style found!")

                # dihedrals
                print "Writing dihedrals..."
                for j, dihedral in enumerate(mol_type.dihedralForceSet.itervalues()):
                    atom1 = mol_type.moleculeSet[0]._atoms[dihedral.atom1 - 1]
                    atomtype1 = atom1.bondtype
                    atom2 = mol_type.moleculeSet[0]._atoms[dihedral.atom2 - 1]
                    atomtype2 = atom2.bondtype
                    atom3 = mol_type.moleculeSet[0]._atoms[dihedral.atom3 - 1]
                    atomtype3 = atom3.bondtype
                    atom4 = mol_type.moleculeSet[0]._atoms[dihedral.atom4 - 1]
                    atomtype4 = atom4.bondtype

                    if isinstance(dihedral, FourierDihedral):
                        style = 'opls'
                        temp = FourierDihedralType(atomtype1, atomtype2,
                                atomtype3, atomtype4, 
                                dihedral.c1, dihedral.c1,
                                dihedral.c2, dihedral.c3)
                        if temp not in dihedral_type_dict:
                            dihedral_type_dict[temp] = dih_type_i
                            dihedral_coeffs.append('{0:d} {1} {2:18.8f} {3:18.8f} {4:18.8f} {5:18.8f}\n'.format(
                                    dih_type_i,
                                    style,
                                    dihedral.c1.in_units_of(self.ENERGY)._value,
                                    dihedral.c2.in_units_of(self.ENERGY)._value,
                                    dihedral.c3.in_units_of(self.ENERGY)._value,
                                    dihedral.c4.in_units_of(self.ENERGY)._value))
                            dih_type_i += 1
                    else:
                        warn("Found unimplemented dihedral type for LAMMPS!")
                        continue

                    dihedral_list.append('{0:-6d} {1:6d} {2:6d} {3:6d} {4:6d} {5:6d}\n'.format(
                            i + j + 1,
                            dihedral_type_dict[temp],
                            dihedral.atom1 + offset,
                            dihedral.atom2 + offset,
                            dihedral.atom3 + offset,
                            dihedral.atom4 + offset))
                    dihedral_style.add(style)
                    if len(dihedral_style) > 1:
                        warn("More than one dihedral style found!")

            # atom specific information
        """
        # Write the actual data file.
        print "    Writing to '{0}'...".format(data_file)
        start = time.time()
        with open(data_file, 'w') as f:
            # front matter
            f.write('Generated by mBuild\n')
            f.write('\n')

            n_atoms = len(atom_list) - 3
            n_bonds = len(bond_list) - 3
            n_angles = len(angle_list) - 3
            n_dihedrals = len(dihedral_list) - 3
            n_impropers = len(improper_list) - 3

            n_atom_types = len(pair_coeff_list) - 3
            n_bond_types = len(bond_coeffs) - 3
            n_angle_types = len(angle_coeffs) - 3
            n_dihedral_types = len(dihedral_coeffs) - 3
            n_improper_types = len(improper_coeffs) - 3

            f.write('{0} atoms\n'.format(n_atoms))
            f.write('{0} bonds\n'.format(n_bonds))
            f.write('{0} angles\n'.format(n_angles))
            f.write('{0} dihedrals\n'.format(n_dihedrals))
            f.write('{0} impropers\n'.format(n_impropers))
            f.write('\n')

            f.write('{0} atom types\n'.format(n_atom_types))
            if n_bond_types > 0:
                f.write('{0} bond types\n'.format(n_bond_types))
            if n_angle_types > 0:
                f.write('{0} angle types\n'.format(n_angle_types))
            if n_dihedral_types > 0:
                f.write('{0} dihedral types\n'.format(n_dihedral_types))
            if n_improper_types > 0:
                f.write('{0} improper types\n'.format(n_improper_types))
            f.write('\n')

            # shifting of box dimensions
            f.write('{0:10.6f} {1:10.6f} xlo xhi\n'.format(x_min, x_max))
            f.write('{0:10.6f} {1:10.6f} ylo yhi\n'.format(y_min, y_max))
            f.write('{0:10.6f} {1:10.6f} zlo zhi\n'.format(z_min, z_max))

            # masses
            for mass in mass_list:
                f.write(mass)

            # forcefield coefficients
            if len(pair_coeff_list) > 3:
                for pair in pair_coeff_list:
                    f.write(pair)
            if len(bond_coeffs) > 3:
                for bond in bond_coeffs:
                    f.write(bond)
            if len(angle_coeffs) > 3:
                for angle in angle_coeffs:
                    f.write(angle)
            if len(dihedral_coeffs) > 3:
                for dihedral in dihedral_coeffs:
                    f.write(dihedral)
            if len(improper_coeffs) > 3:
                for improper in improper_coeffs:
                    f.write(improper)

            # atoms and velocities
            for atom in atom_list:
                f.write(atom)

            # topology
            if len(bond_list) > 3:
                for bond in bond_list:
                    f.write(bond)
            if len(angle_list) > 3:
                for angle in angle_list:
                    f.write(angle)
            if len(dihedral_list) > 3:
                for dihedral in dihedral_list:
                    f.write(dihedral)
            if len(improper_list) > 3:
                for improper in improper_list:
                    f.write(improper)
        print "    Done. ({0:.2f} s)".format(time.time() - start)
        """
        # Write the corresponding input file.
        basename = os.path.splitext(data_file)[0]
        input_filename = '{0}.input'.format(basename)
        with open(input_filename, 'w') as f:
            f.write('units {0}\n'.format(unit_set))
            f.write('atom_style full\n')  # TODO
            f.write('\n')

            f.write('dimension 3\n')  # TODO
            f.write('boundary p p p\n')  # TODO
            f.write('\n')

            # non-bonded
            f.write('pair_style lj/cut/coul/long 9.0 10.0\n')  # TODO: match mdp
            f.write('pair_modify mix geometric\n')  # TODO: match defaults
            f.write('kspace_style pppm 1.0e-5\n')  # TODO: match mdp
            f.write('\n')

            # bonded
            if len(bond_coeffs) > 3:
                f.write('bond_style hybrid {0}\n'.format(
                        " ".join(bond_style)))
            if len(angle_coeffs) > 3:
                f.write('angle_style hybrid {0}\n'.format(
                        " ".join(angle_style)))
            if len(dihedral_coeffs) > 3:
                f.write('dihedral_style hybrid {0}\n'.format(
                        " ".join(dihedral_style)))
            if len(improper_coeffs) > 3:
                f.write('improper_style hybrid {0}\n'.format(
                        " ".join(improper_style)))

            f.write('special_bonds lj {0} {1} {2} coul {3} {4} {5}\n'.format(
                    0.0,
                    0.0,
                    System._sys._ljCorrection,
                    0.0,
                    0.0,
                    System._sys._coulombCorrection))
            f.write('\n')

            # read data
            f.write('read_data {0}\n'.format(os.path.basename(data_file)))
            f.write('\n')

            # output energies
            energy_terms = " ".join(['ebond',
                                     'eangle',
                                     'edihed',
                                     'eimp',
                                     'epair',
                                     'evdwl',
                                     'ecoul',
                                     'elong',
                                     'etail',
                                     'pe'])

            f.write('thermo_style custom {0}\n'.format(energy_terms))
            f.write('\n')

            f.write('run 0\n')
        """
if __name__ == "__main__":
    m = Lammps("bmim.lmp", cwd="..\\lammps_testing")
    print [a for a in m.atoms()]

    print m.boundingbox()

    TreeView(m, bonds=False, angles=True).show()
