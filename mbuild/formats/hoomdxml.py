from __future__ import division
from math import radians
from collections import namedtuple

import numpy as np
import operator

from mbuild.utils.conversion import RB_to_OPLS
from mbuild.utils.geometry import coord_shift
from mbuild.utils.decorators import breaking_change

__all__ = ['write_hoomdxml']


@breaking_change("See PR#463 on github")
def write_hoomdxml(structure, filename, ref_distance=1.0, ref_mass=1.0,
                   ref_energy=1.0, rigid_bodies=None, shift_coords=True,
                   auto_scale=False):
    """Output a HOOMD XML file.

    Parameters
    ----------
    structure : parmed.Structure
        ParmEd structure object
    filename : str
        Path of the output file.
    ref_distance : float, optional, default=1.0, units=nanometers
        Reference distance for conversion to reduced units
    ref_mass : float, optional, default=1.0, units=amu
        Reference mass for conversion to reduced units
    ref_energy : float, optional, default=1.0, units=kJ/mol
        Reference energy for conversion to reduced units
    rigid_bodies : list
        List of rigid body information. An integer value is required
        for each particle corresponding to the number of the rigid body with
        which the particle should be included. A value of None indicates the
        particle is not part of any rigid body.
    shift_coords : bool, optional, default=True
        Shift coordinates from (0, L) to (-L/2, L/2) if necessary.
    auto_scale : bool, optional, default=False
        Automatically use largest sigma value as ref_distance, largest mass value
        as ref_mass and largest epsilon value as ref_energy.

    Returns
    -------
    ReferenceValues : namedtuple
        Values used in scaling

    Example
    -------
    ref_values = ethane.save(filename='ethane-opls.hoomdxml', forcefield_name='oplsaa', auto_scale=True)
    print(ref_values.mass, ref_values.distance, ref_values.energy)


    Notes
    -----
    The following elements are always written:

    * **position** : particle positions
    * **type** : particle types
    * **mass** : particle masses (default 1.0)
    * **charge** : particle charges

    The following elements may be written if applicable:

    * **pair_coeffs** : Pair coefficients for each particle type (assumes a 12-6 LJ pair style). The following information is written for each particle type:

                        * type : particle type
                        * epsilon : LJ epsilon
                        * sigma : LJ sigma

    * **bond_coeffs** : Coefficients for each bond type (assumes a harmonic bond style). The following information is written for each bond type:

                        * type : bond type
                        * k : force constant (units of energy/distance^2)
                        * r0 : bond rest length (units of distance)

    * **bond** : system bonds
    * **angle_coeffs** : Coefficients for each angle type (assumes a harmonic angle style). The following information is written for each angle type:

                         * type : angle type
                         * k : force constant (units of energy/radians^2)
                         * theta : rest angle (units of radians)

    * **angle** : system angles
    * **dihedral_coeffs** : Coefficients for each dihedral type (assumes an OPLS dihedral style). The following information is written for each dihedral type:

                            * type : dihedral type
                            * k1, k2, k3, k4 : force coefficients (units of energy)

    * **dihedral** : system dihedrals

    * **body** : ID of the rigid body to which each particle belongs

    """
    ref_distance *= 10  # Parmed unit hack
    ref_energy /= 4.184  # Parmed unit hack
    forcefield = True
    if structure[0].type == '':
        forcefield = False
    if auto_scale and forcefield:
        ref_mass = max([atom.mass for atom in structure.atoms])
        pair_coeffs = list(set((atom.type,
                                atom.epsilon,
                                atom.sigma) for atom in structure.atoms))
        ref_energy = max(pair_coeffs, key=operator.itemgetter(1))[1]
        ref_distance = max(pair_coeffs, key=operator.itemgetter(2))[2]

    xyz = np.array([[atom.xx, atom.xy, atom.xz] for atom in structure.atoms])
    if shift_coords:
        xyz = coord_shift(xyz, structure.box[:3])

    with open(filename, 'w') as xml_file:
        xml_file.write('<?xml version="1.2" encoding="UTF-8"?>\n')
        xml_file.write('<hoomd_xml version="1.2">\n')
        xml_file.write('<!-- ref_distance (nm) ref_mass (amu) ref_energy (kJ/mol) -->\n')
        xml_file.write('<!-- {} {} {} -->\n'.format(ref_distance, ref_mass, ref_energy))
        xml_file.write('<configuration time_step="0">\n')
        _write_box_information(xml_file, structure, ref_distance)
        _write_particle_information(xml_file, structure, xyz, forcefield,
                ref_distance, ref_mass, ref_energy)
        _write_bond_information(xml_file, structure, ref_distance, ref_energy)
        _write_angle_information(xml_file, structure, ref_energy)
        _write_dihedral_information(xml_file, structure, ref_energy)
        _write_rigid_information(xml_file, rigid_bodies)
        xml_file.write('</configuration>\n')
        xml_file.write('</hoomd_xml>')

    ReferenceValues = namedtuple("ref_values", ["distance", "mass", "energy"])

    return ReferenceValues(ref_distance, ref_mass, ref_energy)


def _write_particle_information(xml_file, structure, xyz, forcefield,
        ref_distance, ref_mass, ref_energy):
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
    e0 = 5.72956500956023e-4 # e^2-mol/kJ-nm, permittivity of free space
    charge_factor = (4.0 * np.pi * e0 * ref_distance * ref_energy)**0.5
    for charge in charges:
        xml_file.write('{}\n'.format(charge/charge_factor))
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
    xml_file.write('<dihedral>\n')
    for dihedral in structure.rb_torsions:
        t1, t2 = dihedral.atom1.type, dihedral.atom2.type,
        t3, t4 = dihedral.atom3.type, dihedral.atom4.type
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
            dihedral_type, *opls_coeffs))
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

    if not all(body is None for body in rigid_bodies):
        xml_file.write('<body>\n')
        for body in rigid_bodies:
            if body is None:
                body = -1
            xml_file.write('{}\n'.format(int(body)))
        xml_file.write('</body>\n')

def _write_box_information(xml_file, structure, ref_distance):
    """Write box information.

    Parameters
    ----------
    xml_file : file object
        The file object of the hoomdxml file being written
    structure : parmed.Structure
        Parmed structure object
    ref_energy : float, default=1.0
        Reference energy for conversion to reduced units

    """
    if np.allclose(structure.box[3:6], np.array([90, 90, 90])):
        box_str = '<box units="sigma"  Lx="{}" Ly="{}" Lz="{}"/>\n'
        xml_file.write(box_str.format(*structure.box[:3] / ref_distance))
    else:
        a, b, c = structure.box[0:3] / ref_distance
        alpha, beta, gamma = np.radians(structure.box[3:6])

        lx = a
        xy = b * np.cos(gamma)
        xz = c * np.cos(beta)
        ly = np.sqrt(b**2 - xy**2)
        yz = (b*c*np.cos(alpha) - xy*xz) / ly
        lz = np.sqrt(c**2 - xz**2 - yz**2)
        box_str = '<box units="sigma"  Lx="{}" Ly="{}" Lz="{}" xy="{}" xz="{}" yz="{}"/>\n'
        xml_file.write(box_str.format(lx, ly, lz, xy, xz, yz))
