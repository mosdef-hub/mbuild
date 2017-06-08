from __future__ import division

from collections import OrderedDict
from copy import deepcopy
from math import floor
import re

import numpy as np
from oset import oset as OrderedSet

from mbuild import Box
from mbuild.utils.io import import_
from mbuild.utils.sorting import natural_sort

__all__ = ['write_gsd']


def write_gsd(structure, filename, ref_distance=1.0, ref_mass=1.0, 
              ref_energy=1.0, rigid_bodies=None):
    """Output a GSD file (HOOMD v2 default data format).

    Parameters
    ----------
    structure : parmed.Structure
        ParmEd Structure object
    filename : str
        Path of the output file.
    ref_distance : float, optional, default=1.0
        Reference distance for conversion to reduced units
    ref_mass : float, optional, default=1.0
        Reference mass for conversion to reduced units
    ref_energy : float, optional, default=1.0
        Reference energy for conversion to reduced units
    rigid_bodies : list of int, optional, default=None
        List of rigid body information. An integer value is required for
        each atom corresponding to the index of the rigid body the particle
        is to be associated with. A value of None indicates the atom is not
        part of a rigid body.

    Notes
    -----
    Force field parameters are not written to the GSD file and must be included
    manually into a HOOMD input script. Work on a HOOMD plugin is underway to
    read force field parameters from a Foyer XML file.

    """

    import_('gsd')
    import gsd.hoomd

    xyz = np.array([[atom.xx, atom.xy, atom.xz] for atom in structure.atoms])

    gsd_file = gsd.hoomd.Snapshot()

    gsd_file.configuration.step = 0
    gsd_file.configuration.dimensions = 3
    gsd_file.configuration.box = np.hstack((structure.box[:3] / ref_distance,
                                            np.zeros(3)))

    _write_particle_information(gsd_file, structure, xyz, ref_distance,
            ref_mass, ref_energy, rigid_bodies)
    if structure.bonds:
        _write_bond_information(gsd_file, structure)
    if structure.angles:
        _write_angle_information(gsd_file, structure)
    if structure.rb_torsions:
        _write_dihedral_information(gsd_file, structure)

    gsd.hoomd.create(filename, gsd_file)

def _write_particle_information(gsd_file, structure, xyz, ref_distance,
        ref_mass, ref_energy, rigid_bodies):
    """Write out the particle information.

    """

    gsd_file.particles.N = len(structure.atoms)
    gsd_file.particles.position = xyz / ref_distance

    types = [atom.name if atom.type == '' else atom.type 
             for atom in structure.atoms]

    unique_types = list(set(types))
    unique_types.sort(key=natural_sort)
    gsd_file.particles.types = unique_types

    typeids = np.array([unique_types.index(t) for t in types])
    gsd_file.particles.typeid = typeids

    masses = np.array([atom.mass for atom in structure.atoms])
    masses[masses==0] = 1.0
    gsd_file.particles.mass = masses / ref_mass

    charges = np.array([atom.charge for atom in structure.atoms])
    e0 = 2.39725e-4
    '''
    Permittivity of free space = 2.39725e-4 e^2/((kcal/mol)(angstrom)),
    where e is the elementary charge
    '''
    charge_factor = (4.0*np.pi*e0*ref_distance*ref_energy)**0.5
    gsd_file.particles.charge = charges / charge_factor

    if rigid_bodies:
        rigid_bodies = [-1 if body is None else body for body in rigid_bodies]
    gsd_file.particles.body = rigid_bodies

def _write_bond_information(gsd_file, structure):
    """Write the bonds in the system.

    Parameters
    ----------
    gsd_file : 
        The file object of the GSD file being written
    structure : parmed.Structure
        Parmed structure object holding system information

    """

    gsd_file.bonds.N = len(structure.bonds)

    unique_bond_types = set()
    for bond in structure.bonds:
        t1, t2 = bond.atom1.type, bond.atom2.type
        if t1 == '' or t2 == '':
            t1, t2 = bond.atom1.name, bond.atom2.name
        t1, t2 = sorted([t1, t2], key=natural_sort)
        try:
            bond_type = ('-'.join((t1, t2)))
        except AttributeError: # no forcefield applied, bond.type is None
            bond_type = ('-'.join((t1, t2)), 0.0, 0.0)
        unique_bond_types.add(bond_type)
    unique_bond_types = sorted(list(unique_bond_types), key=natural_sort)
    gsd_file.bonds.types = unique_bond_types

    bond_typeids = []
    bond_groups = []
    for bond in structure.bonds:
        t1, t2 = bond.atom1.type, bond.atom2.type
        if t1 == '' or t2 == '':
            t1, t2 = bond.atom1.name, bond.atom2.name
        t1, t2 = sorted([t1, t2], key=natural_sort)
        try:
            bond_type = ('-'.join((t1, t2)))
        except AttributeError: # no forcefield applied, bond.type is None
            bond_type = ('-'.join((t1, t2)), 0.0, 0.0)
        bond_typeids.append(unique_bond_types.index(bond_type))
        bond_groups.append((bond.atom1.idx, bond.atom2.idx))

    gsd_file.bonds.typeid = bond_typeids
    gsd_file.bonds.group = bond_groups

def _write_angle_information(gsd_file, structure):
    """Write the angles in the system.

    Parameters
    ----------
    gsd_file : 
        The file object of the GSD file being written
    structure : parmed.Structure
        Parmed structure object holding system information

    """

    gsd_file.angles.N = len(structure.angles)

    unique_angle_types = set()
    for angle in structure.angles:
        t1, t2, t3 = angle.atom1.type, angle.atom2.type, angle.atom3.type
        t1, t3 = sorted([t1, t3], key=natural_sort)
        angle_type = ('-'.join((t1, t2, t3)))
        unique_angle_types.add(angle_type)
    unique_angle_types = sorted(list(unique_angle_types), key=natural_sort)
    gsd_file.angles.types = unique_angle_types

    angle_typeids = []
    angle_groups = []
    for angle in structure.angles:
        t1, t2, t3 = angle.atom1.type, angle.atom2.type, angle.atom3.type
        t1, t3 = sorted([t1, t3], key=natural_sort)
        angle_type = ('-'.join((t1, t2, t3)))
        angle_typeids.append(unique_angle_types.index(angle_type))
        angle_groups.append((angle.atom1.idx, angle.atom2.idx, 
                             angle.atom3.idx))

    gsd_file.angles.typeid = angle_typeids
    gsd_file.angles.group = angle_groups

def _write_dihedral_information(gsd_file, structure):
    """Write the dihedrals in the system.

    Parameters
    ----------
    gsd_file : 
        The file object of the GSD file being written
    structure : parmed.Structure
        Parmed structure object holding system information

    """

    gsd_file.dihedrals.N = len(structure.rb_torsions)

    unique_dihedral_types = set()
    for dihedral in structure.rb_torsions:
        t1, t2 = dihedral.atom1.type, dihedral.atom2.type
        t3, t4 = dihedral.atom3.type, dihedral.atom4.type
        if [t2, t3] == sorted([t2, t3], key=natural_sort):
            dihedral_type = ('-'.join((t1, t2, t3, t4)))
        else:
            dihedral_type = ('-'.join((t4, t3, t2, t1)))
        unique_dihedral_types.add(dihedral_type)
    unique_dihedral_types = sorted(list(unique_dihedral_types), key=natural_sort)
    gsd_file.dihedrals.types = unique_dihedral_types

    dihedral_typeids = []
    dihedral_groups = []
    for dihedral in structure.rb_torsions:
        t1, t2 = dihedral.atom1.type, dihedral.atom2.type
        t3, t4 = dihedral.atom3.type, dihedral.atom4.type
        if [t2, t3] == sorted([t2, t3], key=natural_sort):
            dihedral_type = ('-'.join((t1, t2, t3, t4)))
        else:
            dihedral_type = ('-'.join((t4, t3, t2, t1)))
        dihedral_typeids.append(unique_dihedral_types.index(dihedral_type))
        dihedral_groups.append((dihedral.atom1.idx, dihedral.atom2.idx,
                                dihedral.atom3.idx, dihedral.atom4.idx))

    gsd_file.dihedrals.typeid = dihedral_typeids
    gsd_file.dihedrals.group = dihedral_groups
