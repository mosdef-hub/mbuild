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
              ref_energy=1.0):
    """Output a GSD file (HOOMD default data format).
    
    Parameters
    ----------
    structure : parmed.Structure
        Parmed Structure object
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
    """

    import_('gsd')
    import gsd.hoomd

    forcefield = True
    if structure[0].type == '':
        forcefield = False

    xyz = np.array([[atom.xx, atom.xy, atom.xz] for atom in structure.atoms])

    # Center box at origin and remap coordinates into box
    box.lengths *= 10.0
    box.maxs *= 10.0
    box.mins *= 10.0
    box_init = deepcopy(box)
    box.mins = np.array([-d/2 for d in box_init.lengths])
    box.maxs = np.array([d/2 for d in box_init.lengths])
    
    shift = [box_init.maxs[i] - max for i, max in enumerate(box.maxs)]
    for i, pos in enumerate(xyz):
        for j, coord in enumerate(pos):
            xyz[i, j] -= shift[j]
            rep = floor((xyz[i, j]-box.mins[j]) / box.lengths[j])
            xyz[i, j] -= (rep * box.lengths[j])

    gsd_file = gsd.hoomd.Snapshot()

    gsd_file.configuration.step = 0
    gsd_file.configuration.dimensions = 3
    gsd_file.configuration.box = np.hstack((box.lengths / ref_distance,
                                            np.zeros(3)))

    gsd_file.particles.N = len(structure.atoms)
    gsd_file.particles.position = xyz / ref_distance

    if forcefield:
        types = [atom.type for atom in structure.atoms]
    else:
        types = [atom.name for atom in structure.atoms]

    unique_types = list(set(types))
    unique_types.sort(key=_natural_sort)

    typeids = np.array([unique_types.index(t) for t in types])

    gsd_file.particles.types = unique_types
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

    bonds = [[bond.atom1.idx, bond.atom2.idx] for bond in structure.bonds]
    if bonds:
        bonds = np.asarray(bonds)
        gsd_file.bonds.N = len(bonds)
        if len(structure.bond_types) == 0:
            bond_types = np.zeros(len(bonds),dtype=int)
            gsd_file.bonds.types = ['0']
        else:
            unique_bond_types = dict(enumerate(OrderedSet([(round(bond.type.k,3),
                                                            round(bond.type.req,3)) for bond in structure.bonds])))
            unique_bond_types = OrderedDict([(y,x) for x,y in unique_bond_types.items()])
            bond_types = [unique_bond_types[(round(bond.type.k,3),
                                             round(bond.type.req,3))] for bond in structure.bonds]
            gsd_file.bonds.types = [str(y) for x,y in unique_bond_types.items()]
        gsd_file.bonds.typeid = bond_types
        gsd_file.bonds.group = bonds

    angles = [[angle.atom1.idx,
               angle.atom2.idx, 
               angle.atom3.idx] for angle in structure.angles]
    if angles:
        angles = np.asarray(angles)
        gsd_file.angles.N = len(angles)
        unique_angle_types = dict(enumerate(OrderedSet([(round(angle.type.k,3),
                                                         round(angle.type.theteq,3)) for angle in structure.angles])))
        unique_angle_types = OrderedDict([(y,x) for x,y in unique_angle_types.items()])
        angle_types = [unique_angle_types[(round(angle.type.k,3),
                                           round(angle.type.theteq,3))] for angle in structure.angles]
        gsd_file.angles.types = [str(y) for x,y in unique_angle_types.items()]
        gsd_file.angles.typeid = angle_types
        gsd_file.angles.group = angles

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


def _atoi(text):
    return int(text) if text.isdigit() else text


def _natural_sort(text):
    return [_atoi(a) for a in re.split(r'(\d+)', text)]
