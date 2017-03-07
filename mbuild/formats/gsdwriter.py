from __future__ import division

from collections import OrderedDict
from copy import deepcopy
from math import floor
import re

import numpy as np
from oset import oset as OrderedSet

from mbuild.utils.io import import_

__all__ = ['write_gsd']


def write_gsd(structure, filename, box, ref_distance=1.0, ref_mass=1.0,
              ref_energy=1.0, rigid_bodies=None, popleft_underscore=True):
    """Output a GSD file (HOOMD v2 default data format).
    
    Parameters
    ----------
    structure : parmed.Structure
        ParmEd Structure object
    filename : str
        Path of the output file.
    box : mb.Box
        Box information
    ref_distance : float, optional, default=1.0
        Reference distance for conversion to reduced units
    ref_mass : float, optional, default=1.0
        Reference mass for conversion to reduced units
    ref_energy : float, optional, default=1.0
        Reference energy for conversion to reduced units
    rigid_bodies : list of int, optional, default=None
        List of rigid body information. An integer value is required for
        each atom corresponding to the index of the rigid body the particle
        is to be associated with. A value of -1 indicates the atom is not
        part of a rigid body.
    popleft_underscore : bool, optional, default=True
        If True (default), remove a leading underscore from the particle
        names. This is useful for non-atomistic (e.g. coarse-grained) systems,
        where `Foyer` may need prepending underscores for non-atomistic
        particles.

    Notes
    -----
    Force field information is not written to the GSD file and must be included
    manually into a HOOMD input script. Work on a HOOMD plugin is underway to
    read force field information from a `Foyer` XML file.

    """

    import_('gsd')
    import gsd.hoomd

    forcefield = True
    if structure[0].type == '':
        forcefield = False

    xyz = np.array([[atom.xx, atom.xy, atom.xz] for atom in structure.atoms])

    gsd_file = gsd.hoomd.Snapshot()

    gsd_file.configuration.step = 0
    gsd_file.configuration.dimensions = 3
    gsd_file.configuration.box = np.hstack((box.lengths / ref_distance,
                                            np.zeros(3)))

    _write_particle_information(gsd_file, structure, xyz, forcefield, ref_distance,
            ref_mass, ref_energy, popleft_underscore, rigid_bodies)
    if structure.bonds:
        _write_bond_information(gsd_file, structure)
    if structure.angles:
        _write_angle_information(gsd_file, structure)
    if structure.dihedrals:
        _write_dihedral_information(gsd_file, structure)

        angles = np.asarray(angles)
        gsd_file.angles.N = len(angles)
        unique_angle_types = dict(enumerate(OrderedSet([(round(angle.type.k,3),
                                                         round(angle.type.theteq,3)) for angle in structure.angles])))
        unique_angle_types = OrderedDict([(y,x) for x,y in unique_angle_types.items()])
        angle_types = [unique_angle_types[(round(angle.type.k,3),
                                           round(angle.type.theteq,3))] for angle in structure.angles]

    dihedrals = [[dihedral.atom1.idx,
                  dihedral.atom2.idx,
                  dihedral.atom3.idx,
                  dihedral.atom4.idx] for dihedral in structure.rb_torsions]
    if dihedrals:
        dihedrals = np.asarray(dihedrals)
        gsd_file.dihedrals.N = len(dihedrals)

        unique_dihedral_types = dict(enumerate(OrderedSet([(round(dihedral.type.c0,3),
                                                    round(dihedral.type.c1,3),
                                                    round(dihedral.type.c2,3),
                                                    round(dihedral.type.c3,3),
                                                    round(dihedral.type.c4,3),
                                                    round(dihedral.type.c5,3),
                                                    round(dihedral.type.scee,1),
                                                    round(dihedral.type.scnb,1)) for dihedral in structure.rb_torsions])))
        unique_dihedral_types = OrderedDict([(y,x) for x,y in unique_dihedral_types.items()])
        dihedral_types = [unique_dihedral_types[(round(dihedral.type.c0,3),
                                                 round(dihedral.type.c1,3),
                                                 round(dihedral.type.c2,3),
                                                 round(dihedral.type.c3,3),
                                                 round(dihedral.type.c4,3),
                                                 round(dihedral.type.c5,3),
                                                 round(dihedral.type.scee,1),
                                                 round(dihedral.type.scnb,1))] for dihedral in structure.rb_torsions]
        gsd_file.dihedrals.types = [str(y) for x,y in unique_dihedral_types.items()]
        gsd_file.dihedrals.typeid = dihedral_types
        gsd_file.dihedrals.group = dihedrals

    gsd.hoomd.create(filename, gsd_file)

def _write_particle_information(gsd_file, structure, xyz, ref_distance,
        ref_mass, ref_energy, popleft_underscore, rigid_bodies):
    """Write out the particle information.

    """

    gsd_file.particles.N = len(structure.atoms)
    gsd_file.particles.position = xyz / ref_distance

    if forcefield:
        types = [atom.type for atom in structure.atoms]
    else:
        types = [atom.name for atom in structure.atoms]
    if popleft_underscore:
        for atom_type in types:
            atom_type = atom_type[1:]

    unique_types = list(set(types))
    unique_types.sort(key=_natural_sort)
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
        t1, t2 = sorted([t1, t2])
        try:
            bond_type = ('-'.join((t1, t2)), bond.type.k, bond.type.req)
        except AttributeError: # no forcefield applied, bond.type is None
            bond_type = ('-'.join((t1, t2)), 0.0, 0.0)
        unique_bond_types.add(bond_type)
    unique_bond_types = sorted(list(unique_bond_types), key=_natural_sort)
    gsd_file.bonds.types = unique_bond_types

    bond_typeids = []
    bond_groups = []
    for bond in structure.bonds:
        t1, t2 = bond.atom1.type, bond.atom2.type
        if t1 == '' or t2 == '':
            t1, t2 = bond.atom1.name, bond.atom2.name
        t1, t2 = sorted([t1, t2])
        try:
            bond_type = ('-'.join((t1, t2)), bond.type.k, bond.type.req)
        except AttributeError: # no forcefield applied, bond.type is None
            bond_type = ('-'.join((t1, t2)), 0.0, 0.0)
        bond_typeids.append(unique_bond_types.index(bond_type))
        bond_groups.append((bond.atom1.index, bond.atom2.index))

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
    for angle in structure.angle:
        t1, t2, t3 = angle.atom1.type, angle.atom2.type, angle.atom3.type
        t1, t3 = sorted([t1, t3])
        angle_type = ('-'.join((t1, t2)), angle.type.k, angle.type.req)

        unique_bond_types.add(bond_type)
    unique_bond_types = sorted(list(unique_bond_types), key=_natural_sort)
    gsd_file.bonds.types = unique_bond_types

    gsd_file.angles.types = [str(y) for x,y in unique_angle_types.items()]
    gsd_file.angles.typeid = angle_types
    gsd_file.angles.group = angles

def _write_dihedral_information(gsd_file, structure):

def _atoi(text):
    return int(text) if text.isdigit() else text


def _natural_sort(text):
    return [_atoi(a) for a in re.split(r'(\d+)', text)]
