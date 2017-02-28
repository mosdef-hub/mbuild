from __future__ import division

__all__ = ['write_hoomdxml']


from copy import deepcopy
from math import floor, radians

import numpy as np


def RB_to_OPLS(c0, c1, c2, c3, c4, c5):
    """Converts Ryckaert-Bellemans type dihedrals to OPLS type.

    Parameters
    ----------
    c0, c1, c2, c3, c4, c5 : Ryckaert-Belleman coefficients (in kcal/mol)

    Returns
    -------
    opls_coeffs : np.array, shape=(4,)
        Array containing the OPLS dihedrals coeffs f1, f2, f3, and f4
        (in kcal/mol)
    """

    f1 = (-1.5 * c3) - (2 * c1)
    f2 = c0 + c1 + c3
    f3 = -0.5 * c3
    f4 = -0.25 * c4

    opls_coeffs = np.array([f1, f2, f3, f4])

    return opls_coeffs


def write_hoomdxml(structure, filename, box, rigid_bodies, ref_distance=1.0, 
                   ref_mass=1.0, ref_energy=1.0, wrap_coordinates=True):
    """Output a HOOMD XML file.

    The following elements are always written:

    Position : atomic positions
    Type : atom types
    Mass : atom masses (default 1.0)
    Charge : atom charges

    The following elements may be written if applicable:
    Pair_Coeffs : Pair coefficients for each atom type, assumes a 12-6
                  LJ pair style. The following information is written:
                  type : atom type
                  epsilon : LJ epsilon
                  sigma : LJ sigma
    Bond_Coeffs : Coefficients for each bond type, assumes a harmonic
                  bond style. The following information is written:
                  type : bond type
                  k : force constant (units of energy/distance^2)
                  r0 : bond rest length (units of distance)
    Bond : system bonds
    Angle_Coeffs : Coefficients for each angle type, assumes a harmonic
                   angle style. The following information is written:
                   type : angle type
                   k : force constant (units of energy/radians^2)
                   theta : rest angle (units of radians)
    Angle : system angles
    Dihedral_Coeffs : Coefficients for each dihedral type, assumes an OPLS
                      dihedral style. The following information is written:
                      type : dihedral type
                      k1, k2, k3, k4 : force coefficients (units of energy)
    Dihedral : system dihedrals
    Body : rigid body to which each atom belongs

    Parameters
    ----------
    structure : parmed.Structure
        Parmed structure object
    filename : str
        Path of the output file.
    box : mb.Box
        Box information
    rigid_bodies : list
        List of rigid body information. An integer value is required
        for each atom corresponding to the number of the rigid body with
        which the atom should be included. A value of None indicates the
        atom is not part of any rigid body.
    ref_distance : float, optional, default=1.0
        Reference distance for conversion to reduced units
    ref_mass : float, optional, default=1.0
        Reference mass for conversion to reduced units
    ref_energy : float, optional, default=1.0
        Reference energy for conversion to reduced units
    wrap_coordinates : bool, optional, default=True
        Wrap coordinates of all atoms into the box

    """

    forcefield = True
    if structure[0].type == '':
        forcefield = False

    xyz = np.array([[atom.xx, atom.xy, atom.xz] for atom in structure.atoms])

    # Center box at origin and remap coordinates into box
    box.maxs *= 10.
    box.mins *= 10.
    box_init = deepcopy(box)
    box.mins = np.array([-d/2 for d in box_init.lengths])
    box.maxs = np.array([d/2 for d in box_init.lengths])

    if wrap_coordinates:
        shift = [box_init.maxs[i] - max for i, max in enumerate(box.maxs)]
        for i, pos in enumerate(xyz):
            for j, coord in enumerate(pos):
                xyz[i, j] -= shift[j]
                rep = floor((xyz[i, j] - box.mins[j]) / box.lengths[j])
                xyz[i, j] -= (rep * box.lengths[j])

    with open(filename, 'w') as xml_file:
        xml_file.write('<?xml version="1.2" encoding="UTF-8"?>\n')
        xml_file.write('<hoomd_xml version="1.2">\n')
        xml_file.write('<configuration time_step="0">\n')
        xml_file.write('<box units="sigma"  Lx="{}" Ly="{}" Lz="{}"/>\n'.format(*box.lengths/ref_distance))

        xml_file.write('<position units="sigma" num="{}">\n'.format(len(structure.atoms)))

        for pos in xyz:
            xml_file.write('{}\t{}\t{}\n'.format(*pos/ref_distance))
        xml_file.write('</position>\n')

        if forcefield:
            types = [atom.type for atom in structure.atoms]
        else:
            types = [atom.name for atom in structure.atoms]
        xml_file.write('<type>\n')
        for atom_type in types:
            xml_file.write('{}\n'.format(atom_type))
        xml_file.write('</type>\n')

        masses = [atom.mass for atom in structure.atoms]
        xml_file.write('<mass>\n')
        for mass in masses:
            if mass == 0:
                mass = 1.0
            xml_file.write('{}\n'.format(mass/ref_mass))
        xml_file.write('</mass>\n')

        charges = [atom.charge for atom in structure.atoms]
        xml_file.write('<charge>\n')
        for charge in charges:
            xml_file.write('{}\n'.format(charge))
        xml_file.write('</charge>\n')

        if forcefield:
            pair_coeffs = list(set((atom.type,
                                    atom.epsilon,
                                    atom.sigma) for atom in structure.atoms))
            pair_coeffs.sort(key=lambda pair_type: pair_type[0])
            xml_file.write('<pair_coeffs>\n')
            for param_set in pair_coeffs:
                xml_file.write('{}\t{:.4f}\t{:.4f}\n'.format(param_set[0], param_set[1]/ref_energy,
                                                             param_set[2]/ref_distance))
            xml_file.write('</pair_coeffs>\n')

        bonds = [[bond.atom1.idx, bond.atom2.idx] for bond in structure.bonds]
        if bonds:
            if len(structure.bond_types) == 0:
                bond_types = np.zeros(len(bonds), dtype=int)
            else:
                unique_bond_types = dict(enumerate(set([(round(bond.type.k, 3),
                                                         round(bond.type.req, 3)) for bond in structure.bonds])))
                unique_bond_types = {y: x for x, y in unique_bond_types.items()}
                bond_types = [unique_bond_types[(round(bond.type.k, 3),
                                                 round(bond.type.req, 3))] for bond in structure.bonds]
                xml_file.write('<bond_coeffs>\n')
                for params, idx in unique_bond_types.items():
                    xml_file.write('{}\t{}\t{}\n'.format(idx, ((params[0]*2.) / ref_energy)*(ref_distance**2.), params[1]/ref_distance))
                xml_file.write('</bond_coeffs>\n')
            xml_file.write('<bond>\n')
            for idx, bond in enumerate(bonds):
                xml_file.write('{}\t{}\t{}\n'.format(bond_types[idx], *bond))
            xml_file.write('</bond>\n')

        angles = [[angle.atom1.idx,
                   angle.atom2.idx,
                   angle.atom3.idx] for angle in structure.angles]
        if angles:
            unique_angle_types = dict(enumerate(set([(round(angle.type.k, 3),
                                                      round(angle.type.theteq, 3)) for angle in structure.angles])))
            unique_angle_types = {y: x for x, y in unique_angle_types.items()}
            angle_types = [unique_angle_types[(round(angle.type.k, 3),
                                               round(angle.type.theteq, 3))] for angle in structure.angles]
            xml_file.write('<angle_coeffs>\n')
            for params, idx in unique_angle_types.items():
                xml_file.write('{}\t{}\t{:.5f}\n'.format(idx, (params[0]*2.)/ref_energy,
                                                         radians(params[1])))
            xml_file.write('</angle_coeffs>\n')
            xml_file.write('<angle>\n')
            for idx, angle in enumerate(angles):
                xml_file.write('{}\t{}\t{}\t{}\n'.format(angle_types[idx], *angle))
            xml_file.write('</angle>\n')

        dihedrals = [[dihedral.atom1.idx,
                      dihedral.atom2.idx,
                      dihedral.atom3.idx,
                      dihedral.atom4.idx] for dihedral in structure.rb_torsions]
        if dihedrals:
            unique_dihedral_types = dict(enumerate(set([(round(dihedral.type.c0, 3),
                                                         round(dihedral.type.c1, 3),
                                                         round(dihedral.type.c2, 3),
                                                         round(dihedral.type.c3, 3),
                                                         round(dihedral.type.c4, 3),
                                                         round(dihedral.type.c5, 3),
                                                         round(dihedral.type.scee, 1),
                                                         round(dihedral.type.scnb, 1)) for dihedral in structure.rb_torsions])))
            unique_dihedral_types = {y: x for x, y in unique_dihedral_types.items()}
            dihedral_types = [unique_dihedral_types[(round(dihedral.type.c0, 3),
                                                     round(dihedral.type.c1, 3),
                                                     round(dihedral.type.c2, 3),
                                                     round(dihedral.type.c3, 3),
                                                     round(dihedral.type.c4, 3),
                                                     round(dihedral.type.c5, 3),
                                                     round(dihedral.type.scee, 1),
                                                     round(dihedral.type.scnb, 1))] for dihedral in structure.rb_torsions]
            xml_file.write('<dihedral_coeffs>\n')
            for params, idx in unique_dihedral_types.items():
                opls_coeffs = RB_to_OPLS(params[0],
                                         params[1],
                                         params[2],
                                         params[3],
                                         params[4],
                                         params[5])
                opls_coeffs /= ref_energy
                xml_file.write('{}\t{:.5f}\t{:.5f}\t{:.5f}\t{:.5f}\n'.format(idx, *opls_coeffs))
            xml_file.write('</dihedral_coeffs>\n')
            xml_file.write('<dihedral>\n')
            for idx, dihedral in enumerate(dihedrals):
                xml_file.write('{}\t{}\t{}\t{}\t{}\n'.format(dihedral_types[idx],
                                                             *dihedral))
            xml_file.write('</dihedral>\n')

        if not all(body is None for body in rigid_bodies):
            xml_file.write('<body>\n')
            for body in rigid_bodies:
                if body is None:
                    body = -1
                xml_file.write('{}\n'.format(int(body)))
            xml_file.write('</body>\n')

        xml_file.write('</configuration>\n')
        xml_file.write('</hoomd_xml>')
