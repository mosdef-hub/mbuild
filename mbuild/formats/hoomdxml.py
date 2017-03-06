from __future__ import division
from copy import deepcopy
from math import floor, radians

import numpy as np

__all__ = ['write_hoomdxml']


def write_hoomdxml(structure, filename, box, ref_distance=1.0, ref_mass=1.0,
                   ref_energy=1.0, rigid_bodies=None, wrap_coordinates=True,
                   popleft_underscore=True):
    """Output a HOOMD XML file.

    Parameters
    ----------
    structure : parmed.Structure
        Parmed structure object
    filename : str
        Path of the output file.
    box : mb.Box
        Box information
    ref_distance : float, default=1.0
        Reference distance for conversion to reduced units
    ref_mass : float, default=1.0
        Reference mass for conversion to reduced units
    ref_energy : float, default=1.0
        Reference energy for conversion to reduced units
    rigid_bodies : list, default=None
        List of rigid body information following the HOOMD XML format.
        An integer value is required for each atom corresponding to the
        number of the rigid body with which the atom should be included.
        A value of -1 indicates the atom is not part of a rigid body.
    popleft_underscore : bool, default=True
        If True (default), remove a leading underscore from the particle names.
        This is useful for non-atomistic (e.g., coarse-grained) systems, where
        `Foyer` may need prepending underscores for non-atomistic particles.

    Notes
    -----
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

    """

    forcefield = True
    if structure[0].type == '':
        forcefield = False
    xyz = np.array([[atom.xx, atom.xy, atom.xz] for atom in structure.atoms])

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
        xml_file.write(
                '<box units="sigma"  Lx="{}" Ly="{}" Lz="{}"/>\n'.format(
                    *box.lengths/ref_distance))
        _write_particle_information(xml_file, structure, xyz, forcefield,
                ref_distance, ref_mass, ref_energy, popleft_underscore)
        _write_bond_information(xml_file, structure, ref_distance, ref_energy)
        _write_angle_information(xml_file, structure, ref_energy)
        _write_dihedral_information(xml_file, structure, ref_energy)
        _write_rigid_information(xml_file, rigid_bodies)
        xml_file.write('</configuration>\n')
        xml_file.write('</hoomd_xml>')


def _write_particle_information(xml_file, structure, xyz, forcefield,
        ref_distance, ref_mass, ref_energy, popleft_underscore):
    """Write out the particle information.

    Parameters
    ----------
    xml_file : file object
        The file object of the hoomdxml file being written
    structure : parmed.Structure
        Parmed structure object
    xyz : np.ndarray, shape=(n,3), dtype=float
        The particle positions to be written.
    forcefield : bool
        If True, write the particle "type". Write the particle "name" otherwise.
    ref_distance : float, default=1.0
        Reference distance for conversion to reduced units
    ref_mass : float, default=1.0
        Reference mass for conversion to reduced units
    ref_energy : float, default=1.0
        Reference energy for conversion to reduced units
    popleft_underscore : bool, default=True
        If True, remove a leading underscore from the particle names.
        This is useful for non-atomistic (e.g., coarse-grained) systems, where
        `Foyer` may need prepending underscores for non-atomistic particles.

    """

    xml_file.write('<position units="sigma" num="{}">\n'.format(xyz.shape[0]))
    for pos in xyz:
        xml_file.write('{}\t{}\t{}\n'.format(*pos/ref_distance))
    xml_file.write('</position>\n')
    if forcefield:
        types = [atom.type for atom in structure.atoms]
    else:
        types = [atom.name for atom in structure.atoms]

    xml_file.write('<type>\n')
    for atom_type in types:
        if popleft_underscore and atom_type.startswith('_'):
            atom_type = atom_type[1:]
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
            xml_file.write('{}\t{:.4f}\t{:.4f}\n'.format(
                param_set[0], param_set[1]/ref_energy,
                param_set[2]/ref_distance))
        xml_file.write('</pair_coeffs>\n')


def _write_bond_information(xml_file, structure, ref_distance, ref_energy):
    """Write the bonds in the system.

    Parameters
    ----------
    xml_file : file object
        The file object of the hoomdxml file being written
    structure : parmed.Structure
        Parmed structure object
    ref_distance : float, default=1.0
        Reference distance for conversion to reduced units
    ref_energy : float, default=1.0
        Reference energy for conversion to reduced units

    """

    unique_bond_types = set()
    xml_file.write('<bond>\n')
    for bond in structure.bonds:
        t1, t2 = bond.atom1.type, bond.atom2.type
        if t1 == '' or t2 == '':
            t1, t2 = bond.atom1.name, bond.atom2.name
        t1, t2 = sorted([t1, t2])
        try:
            bond_type = ('-'.join((t1, t2)), bond.type.k, bond.type.req)
        except AttributeError:  # no forcefield applied, bond.type is None
            bond_type = ('-'.join((t1, t2)), 0.0, 0.0)
        unique_bond_types.add(bond_type)
        xml_file.write('{} {} {}\n'.format(
            bond_type[0], bond.atom1.idx, bond.atom2.idx))
    xml_file.write('</bond>\n')
    xml_file.write('<bond_coeffs>\n')
    xml_file.write('<!-- type k r_eq -->\n')
    for bond_type, k, req in unique_bond_types:
        xml_file.write('{} {} {}\n'.format(bond_type,
            k * 2.0 / ref_energy * ref_distance**2.0, req/ref_distance))
    xml_file.write('</bond_coeffs>\n')


def _write_angle_information(xml_file, structure, ref_energy):
    """Write the angles in the system.

    Parameters
    ----------
    xml_file : file object
        The file object of the hoomdxml file being written
    structure : parmed.Structure
        Parmed structure object
    ref_energy : float, default=1.0
        Reference energy for conversion to reduced units

    """

    unique_angle_types = set()
    xml_file.write('<angle>\n')
    for angle in structure.angles:
        t1, t2, t3 = angle.atom1.type, angle.atom2.type, angle.atom3.type
        t1, t3 = sorted([t1, t3])
        angle_type = ('-'.join((t1, t2, t3)), angle.type.k, angle.type.theteq)
        unique_angle_types.add(angle_type)
        xml_file.write('{} {} {} {}\n'.format(
            angle_type[0], angle.atom1.idx, angle.atom2.idx, angle.atom3.idx))
    xml_file.write('</angle>\n')
    xml_file.write('<angle_coeffs>\n')
    xml_file.write('<!-- type k theta_eq -->\n')
    for angle_type, k, teq in unique_angle_types:
        xml_file.write('{} {} {}\n'.format(angle_type,
            k * 2.0 / ref_energy, radians(teq)))
    xml_file.write('</angle_coeffs>\n')


def _write_dihedral_information(xml_file, structure, ref_energy):
    """Write dihedrals in the system.

    Parameters
    ----------
    xml_file : file object
        The file object of the hoomdxml file being written
    structure : parmed.Structure
        Parmed structure object
    ref_energy : float, default=1.0
        Reference energy for conversion to reduced units

    """

    unique_dihedral_types = set()
    xml_file.write('<dihderal>\n')
    for dihedral in structure.rb_torsions:
        t1, t2 = dihedral.atom1.type, dihedral.atom2.type,
        t3, t4 = dihedral.atom3.type, dihedral.atom4.type
        st2, st3 = sorted([t2, t3])  # sort based on the middle 2
        if [t2, t3] == sorted([t2, t3]):
            types_in_dihedral = '-'.join((t1, t2, t3, t4))
        else:
            types_in_dihedral = '-'.join((t4, t3, t2, t1))
        dihedral_type = (types_in_dihedral, dihedral.type.c0,
        dihedral.type.c1, dihedral.type.c2, dihedral.type.c3, dihedral.type.c4,
        dihedral.type.c5, dihedral.type.scee, dihedral.type.scnb)
        unique_dihedral_types.add(dihedral_type)
        xml_file.write('{} {} {} {} {}\n'.format(
            dihedral_type[0], dihedral.atom1.idx, dihedral.atom2.idx,
            dihedral.atom3.idx, dihedral.atom4.idx))
    xml_file.write('</dihedral>\n')
    xml_file.write('<dihedral_coeffs>\n')
    xml_file.write('<!-- type k1 k2 k3 k4 -->\n')
    for dihedral_type, c0, c1, c2, c3, c4, c5, scee, scnb in unique_dihedral_types:
        opls_coeffs = RB_to_OPLS(c0, c1, c2, c3, c4, c5)
        opls_coeffs /= ref_energy
        xml_file.write('{} {:.5f} {:.5f} {:.5f} {:.5f}\n'.format(
            dihedral_type[0], *opls_coeffs))
    xml_file.write('</dihedral_coeffs>\n')


def _write_rigid_information(xml_file, rigid_bodies):
    """Write rigid body information.

    Parameters
    ----------
    xml_file : file object
        The file object of the hoomdxml file being written
    rigid_bodies : list, len=n_particles
        The rigid body that each particle belongs to (-1 for none)

    """

    if rigid_bodies is not None:
        xml_file.write('<body>\n')
        for body in rigid_bodies:
            xml_file.write('{}\n'.format(int(body)))
        xml_file.write('</body>\n')


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
    return np.array([f1, f2, f3, f4])
