from __future__ import division

__all__ = ['write_gsd']


from copy import deepcopy
from math import floor
import numpy as np
import gsd.hoomd
from .hoomdxml import RB_to_OPLS
from .ff_to_json import write_forcefield

from oset import oset as OrderedSet
from collections import OrderedDict


def write_gsd(structure, filename, forcefield, box, ref_distance=1.0, ref_mass=1.0,
              ref_energy=1.0, write_ff=True):
    """Output a GSD file (HOOMD default data format).
    
    Parameters
    ----------
    structure : parmed.GromacsTopologyFile
        Parmed structure object
    filename : str
        Path of the output file.
    box : mb.Box
        Box information to save to XML file
    forcefield : str, default=None
        Name of the force field to be applied to the compound
    ref_distance : float, default=1.0
        Reference distance for conversion to reduced units
    ref_mass : float, default=1.0
        Reference mass for conversion to reduced units
    ref_energy : float, default=1.0
        Reference energy for conversion to reduced units
    write_ff : boolean, default=True
        Write forcefield parameters to a JSON file, 'parameters.json'
    """

    xyz = np.array([[atom.xx,atom.xy,atom.xz] for atom in structure.atoms])

    # Center box at origin and remap coordinates into box
    box.lengths *= 10.0
    box.maxs *= 10.0
    box.mins *= 10.0
    box_init = deepcopy(box)
    box.mins = np.array([-d/2 for d in box_init.lengths])
    box.maxs = np.array([d/2 for d in box_init.lengths])
    
    shift = [box_init.maxs[i] - max for i,max in enumerate(box.maxs)]
    for i,pos in enumerate(xyz):
        for j,coord in enumerate(pos):
            xyz[i,j] -= shift[j]
            rep = floor((xyz[i,j]-box.mins[j]) / box.lengths[j])
            xyz[i,j] -= (rep * box.lengths[j])

    gsd_file = gsd.hoomd.Snapshot()

    gsd_file.configuration.step = 0
    gsd_file.configuration.dimensions = 3
    gsd_file.configuration.box = np.hstack((box.lengths / ref_distance,np.zeros(3)))

    gsd_file.particles.N = len(structure.atoms)
    gsd_file.particles.position = xyz / ref_distance

    if forcefield:
        types = [atom.type for atom in structure.atoms]
    else:
        types = [atom.name for atom in structure.atoms]
    ntypes = 0
    mapping = OrderedDict()
    for t in types:
        if t not in mapping:
            mapping[t] = ntypes
            ntypes = ntypes+1
    typeids = np.array([mapping[t] for t in types])

    unique_types = OrderedSet(types)
    for unique_type in unique_types:
        ref_atom = structure.atoms[types.index(unique_type)]

    gsd_file.particles.types = types
    gsd_file.particles.typeids = typeids
    
    masses = np.array([atom.mass for atom in structure.atoms])
    masses[masses==0] = 1.0
    gsd_file.particles.mass = masses / ref_mass

    charges = np.arrage([atom.charge for atom in structure.atoms])
    charge_factor = (4.0*np.pi*e0*ref_distance*ref_energy)**0.5
    gsd_file.particles.charge = charges / charge_factor

    bonds = [[bond.atom1.idx, bond.atom2.idx] for bond in structure.bonds]
    if bonds:
        bonds = np.asarray(bonds)
        gsd_file.bonds.N = len(bonds)

        if len(structure.bond_types) == 0:
            bond_types = np.zeros(len(bonds),dtype=int)
        else:
            all_bond_types = OrderedDict(enumerate(set([(round(bond.type.k,3),
                                                         round(bond.type.req,3)) for bond in structure.bonds])))
            all_bond_types = OrderedDict([(y,x) for x,y in all_bond_types.items()])
            bond_types = np.array([all_bond_types[(round(bond.type.k,3),round(bond.type.req,3))] for bond in structure.bonds])
        gsd_file.bonds.typeid = bond_types
        gsd_file.bonds.types = [str(t) for t in bond_types]
        gsd_file.bonds.group = bonds

    angles = [[angle.atom1.idx,
               angle.atom2.idx, 
               angle.atom3.idx] for angle in structure.angles]
    if angles:
        angles = np.asarray(angles)
        gsd_file.angles.N = len(angles)
        all_angle_types = OrderedDict(enumerate(set([(round(angle.type.k,3), 
                                                      round(angle.type.theteq,3)) for angle in structure.angles])))
        all_angle_types = OrderedDict([(y,x) for x,y in all_angle_types.items()])
        angle_types = np.array([all_angle_types[(round(angle.type.k,3),round(angle.type.theteq,3))] for angle in structure.angles])
        gsd_file.angles.typeid = angle_types
        gsd_file.angles.types = [str(t) for t in angle_types]
        gsd_file.angles.group = angles

    dihedrals = [[dihedral.atom1.idx,
                  dihedral.atom2.idx,
                  dihedral.atom3.idx,
                  dihedral.atom4.idx] for dihedral in structure.rb_torsions]
    if dihedrals:
        dihedrals = np.asarray(dihedrals)
        gsd_file.dihedrals.N = len(dihedrals)
        all_dihedral_types = OrderedDict(enumerate(set([(round(dihedral.type.c0,3),
                                                         round(dihedral.type.c1,3),
                                                         round(dihedral.type.c2,3),
                                                         round(dihedral.type.c3,3),
                                                         round(dihedral.type.c4,3),
                                                         round(dihedral.type.c5,3),
                                                         round(dihedral.type.scee,1),
                                                         round(dihedral.type.scnb,1)) for dihedral in structure.rb_torsions])))
        all_dihedral_types = OrderedDict([(y,x) for x,y in all_dihedral_types.items()])
        dihedral_types = np.array([all_dihedral_types[(round(dihedral.type.c0,3),
                                                round(dihedral.type.c1,3),
                                                round(dihedral.type.c2,3),
                                                round(dihedral.type.c3,3),
                                                round(dihedral.type.c4,3),
                                                round(dihedral.type.c5,3),
                                                round(dihedral.type.scee,1),
                                                round(dihedral.type.scnb,1))] for dihedral in structure.rb_torsions])
        gsd_file.dihedrals.typeid = dihedral_types
        gsd_file.dihedrals.types = [str(t) for t in dihedral_types]
        gsd_file.dihedrals.group = dihedrals

    gsd.hoomd.create(filename, gsd_file)

    if write_ff:
        write_forcefield(structure, 'parameters.json', ref_distance=ref_distance, ref_energy=ref_energy)
