from collections import namedtuple
import numpy as np
import parmed as pmd

import mbuild as mb
from mbuild.utils.sorting import natural_sort
from mbuild.utils.geometry import coord_shift
from mbuild.utils.io import import_

hoomd = import_("hoomd")
hoomd.data = import_("hoomd.data")

__all__ = ['to_hoomdsnapshot']

def to_hoomdsnapshot(structure,  ref_distance=1.0, ref_mass=1.0,
              ref_energy=1.0, rigid_bodies=None, shift_coords=True,
              write_special_pairs=True, auto_scale=False, parmed_kwargs={}):
    """Convert mb.Compound or parmed.Structure to hoomd.data.Snapshot

    Parameters
    ----------
    structure : parmed.Structure
        ParmEd Structure object
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
    shift_coords : bool, optional, default=True
        Shift coordinates from (0, L) to (-L/2, L/2) if necessary.
    auto_scale : bool, optional, default=False
        Automatically use largest sigma value as ref_distance,
        largest mass value as ref_mass
        and largest epsilon value as ref_energy
    write_special_pairs : bool, optional, default=True
        Writes out special pair information necessary to correctly use 
        the OPLS fudged 1,4 interactions in HOOMD.

    Returns
    -------
    hoomd_snapshot : hoomd.data.Snapshot
    ReferenceValues : namedtuple
        Values used in scaling


    Notes
    -----
    Force field parameters are not written to the hoomd_snapshot 

    """
    if not isinstance(structure, (mb.Compound, pmd.Structure)):
        raise ValueError("You are trying to create a hoomd.Snapshot from " +
                "{} ".format(type(structure)) + 
                "please pass mb.Compound or pmd.Structure")
    elif isinstance(structure, mb.Compound):
        structure = structure.to_parmed(**parmed_kwargs)


    if not hoomd.context.current:
        hoomd.context.initialize("")

    if auto_scale:
        ref_mass = max([atom.mass for atom in structure.atoms])
        pair_coeffs = list(set((atom.type,
                                atom.epsilon,
                                atom.sigma) for atom in structure.atoms))
        ref_energy = max(pair_coeffs, key=operator.itemgetter(1))[1]
        ref_distance = max(pair_coeffs, key=operator.itemgetter(2))[2]

    ReferenceValues = namedtuple("ref_values", ["distance", "mass", "energy"])
    ref_values = ReferenceValues(ref_distance, ref_mass, ref_energy)

    xyz = np.array([[atom.xx, atom.xy, atom.xz] for atom in structure.atoms])
    if shift_coords:
        xyz = coord_shift(xyz, structure.box[:3])


    # Get box information
    if np.allclose(structure.box[3:6], np.array([90, 90, 90])):
        lx,ly,lz = structure.box[:3]/ref_distance
        xy,xz,yz = 0,0,0
    else:
        a, b, c = structure.box[0:3] / ref_distance
        alpha, beta, gamma = np.radians(structure.box[3:6])

        lx = a
        xy = b * np.cos(gamma)
        xz = c * np.cos(beta)
        ly = np.sqrt(b**2 - xy**2)
        yz = (b*c*np.cos(alpha) - xy*xz) / ly
        lz = np.sqrt(c**2 - xz**2 - yz**2)

    n_particles, scaled_positions, unique_types, typeids, scaled_mass, scaled_charges, rigid_bodies = _parse_particle_information(structure, xyz, ref_distance, ref_mass, ref_energy, rigid_bodies)
    n_bonds, unique_bond_types, bond_typeids, bond_groups = _parse_bond_information(structure)
    n_angles, unique_angle_types, angle_typeids, angle_groups = _parse_angle_information(structure)
    n_dihedrals, unique_dihedral_types, dihedral_typeids, dihedral_groups = _parse_dihedral_information(structure)
    n_impropers, unique_improper_types, improper_typeids, improper_groups = _parse_improper_information(structure)
    pair_types, pair_typeid, pairs, n_pairs = _parse_pair_information(structure)

    hoomd_snapshot = hoomd.data.make_snapshot(N=n_particles,
            box=hoomd.data.boxdim(Lx=lx, Ly=ly, Lz=lz, xy=xy, xz=xz, yz=yz),
            particle_types=unique_types, bond_types=unique_bond_types,
            angle_types=unique_angle_types, dihedral_types=unique_dihedral_types,
            improper_types=unique_improper_types, pair_types=pair_types)

    hoomd_snapshot.particles.resize(n_particles)
    hoomd_snapshot.particles.position[:] = scaled_positions
    hoomd_snapshot.particles.types[:] = unique_types
    hoomd_snapshot.particles.typeid[:] = typeids
    hoomd_snapshot.particles.mass[:] = scaled_mass 
    hoomd_snapshot.particles.charge[:] = scaled_charges
    hoomd_snapshot.particles.body[:] = rigid_bodies

    if n_bonds > 0:
        hoomd_snapshot.bonds.resize(n_bonds)
        hoomd_snapshot.bonds.types[:] = unique_bond_types
        hoomd_snapshot.bonds.typeid[:] = bond_typeids
        hoomd_snapshot.bonds.group[:] = bond_groups

    if n_angles > 0:
        hoomd_snapshot.angles.resize(n_angles)
        hoomd_snapshot.angles.types[:] = unique_angle_types
        hoomd_snapshot.angles.typeid[:] = angle_typeids
        hoomd_snapshot.angles.group[:] = np.reshape(angle_groups, (-1, 3))

    if n_dihedrals > 0:
        hoomd_snapshot.dihedrals.resize(n_dihedrals)
        hoomd_snapshot.dihedrals.types[:] = unique_dihedral_types
        hoomd_snapshot.dihedrals.typeid[:] = dihedral_typeids
        hoomd_snapshot.dihedrals.group[:] = np.reshape(dihedral_groups, (-1,4))

    if n_impropers > 0:
        hoomd_snapshot.impropers.resize(n_impropers)
        hoomd_snapshot.impropers.types[:] = unique_improper_types
        hoomd_snapshot.impropers.typeid[:] = improper_typeids
        hoomd_snapshot.impropers.group[:] = np.reshape(improper_groups, (-1,4))

    if n_pairs > 0:
        hoomd_snapshot.pairs.resize(n_pairs)
        hoomd_snapshot.pairs.types[:] = pair_types
        hoomd_snapshot.pairs.typeid[:] = pair_typeid
        hoomd_snapshot.pairs.group[:] = np.reshape(pairs, (-1,2))

    return hoomd_snapshot, ref_values

def _parse_particle_information(structure, xyz, ref_distance,
        ref_mass, ref_energy, rigid_bodies):
    """Write out the particle information.

    """

    n_particles = len(structure.atoms)
    scaled_positions = xyz/ref_distance

    types = [atom.name if atom.type == '' else atom.type
             for atom in structure.atoms]

    unique_types = list(set(types))
    unique_types.sort(key=natural_sort)

    typeids = np.array([unique_types.index(t) for t in types])

    masses = np.array([atom.mass for atom in structure.atoms])
    masses[masses==0] = 1.0
    scaled_mass = masses / ref_mass

    charges = np.array([atom.charge for atom in structure.atoms])
    e0 = 2.39725e-4
    '''
    Permittivity of free space = 2.39725e-4 e^2/((kcal/mol)(angstrom)),
    where e is the elementary charge
    '''
    charge_factor = (4.0*np.pi*e0*ref_distance*ref_energy)**0.5
    scaled_charges = charges / charge_factor

    if rigid_bodies:
        rigid_bodies = [-1 if body is None else body for body in rigid_bodies]
    else:
        rigid_bodies = [-1 for _ in structure.atoms]

    return (n_particles, scaled_positions, unique_types, typeids, 
            scaled_mass, scaled_charges, rigid_bodies)

def _parse_pair_information(structure):
    pair_types = []
    pair_typeid = [] 
    pairs = [] 
    for ai in structure.atoms:
        for aj in ai.dihedral_partners:
            #make sure we don't double add
            if ai.idx > aj.idx:
                ps = '-'.join(sorted([ai.type, aj.type]))
                if ps not in pair_types:
                    pair_types.append(ps)
                pair_typeid.append(pair_types.index(ps))
                pairs.append((ai.idx, aj.idx))
    n_pairs = len(pairs)

    return pair_types, pair_typeid, pairs, n_pairs

def _parse_bond_information(structure):
    n_bonds = len(structure.bonds)
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


    return n_bonds, unique_bond_types, bond_typeids, bond_groups

def _parse_angle_information(structure):
    n_angles = len(structure.angles)

    unique_angle_types = set()
    for angle in structure.angles:
        t1, t2, t3 = angle.atom1.type, angle.atom2.type, angle.atom3.type
        t1, t3 = sorted([t1, t3], key=natural_sort)
        angle_type = ('-'.join((t1, t2, t3)))
        unique_angle_types.add(angle_type)
    unique_angle_types = sorted(list(unique_angle_types), key=natural_sort)

    angle_typeids = []
    angle_groups = []
    for angle in structure.angles:
        t1, t2, t3 = angle.atom1.type, angle.atom2.type, angle.atom3.type
        t1, t3 = sorted([t1, t3], key=natural_sort)
        angle_type = ('-'.join((t1, t2, t3)))
        angle_typeids.append(unique_angle_types.index(angle_type))
        angle_groups.append((angle.atom1.idx, angle.atom2.idx,
                             angle.atom3.idx))

    return n_angles, unique_angle_types, angle_typeids, angle_groups

def _parse_dihedral_information(structure):
    n_dihedrals = len(structure.rb_torsions + structure.dihedrals)

    unique_dihedral_types = set()
    for dihedral in structure.rb_torsions+structure.dihedrals:
        t1, t2 = dihedral.atom1.type, dihedral.atom2.type
        t3, t4 = dihedral.atom3.type, dihedral.atom4.type
        if [t2, t3] == sorted([t2, t3], key=natural_sort):
            dihedral_type = ('-'.join((t1, t2, t3, t4)))
        else:
            dihedral_type = ('-'.join((t4, t3, t2, t1)))
        unique_dihedral_types.add(dihedral_type)
    unique_dihedral_types = sorted(list(unique_dihedral_types), key=natural_sort)

    dihedral_typeids = []
    dihedral_groups = []
    for dihedral in structure.rb_torsions+structure.dihedrals:
        t1, t2 = dihedral.atom1.type, dihedral.atom2.type
        t3, t4 = dihedral.atom3.type, dihedral.atom4.type
        if [t2, t3] == sorted([t2, t3], key=natural_sort):
            dihedral_type = ('-'.join((t1, t2, t3, t4)))
        else:
            dihedral_type = ('-'.join((t4, t3, t2, t1)))
        dihedral_typeids.append(unique_dihedral_types.index(dihedral_type))
        dihedral_groups.append((dihedral.atom1.idx, dihedral.atom2.idx,
                                dihedral.atom3.idx, dihedral.atom4.idx))


    return n_dihedrals, unique_dihedral_types, dihedral_typeids, dihedral_groups

def _parse_improper_information(structure):
    n_dihedrals = len(structure.impropers)

    unique_dihedral_types = set()
    for dihedral in structure.impropers:
        t1, t2 = dihedral.atom1.type, dihedral.atom2.type
        t3, t4 = dihedral.atom3.type, dihedral.atom4.type
        if [t2, t3] == sorted([t2, t3], key=natural_sort):
            dihedral_type = ('-'.join((t1, t2, t3, t4)))
        else:
            dihedral_type = ('-'.join((t4, t3, t2, t1)))
        unique_dihedral_types.add(dihedral_type)
    unique_dihedral_types = sorted(list(unique_dihedral_types), key=natural_sort)

    dihedral_typeids = []
    dihedral_groups = []
    for dihedral in structure.impropers:
        t1, t2 = dihedral.atom1.type, dihedral.atom2.type
        t3, t4 = dihedral.atom3.type, dihedral.atom4.type
        if [t2, t3] == sorted([t2, t3], key=natural_sort):
            dihedral_type = ('-'.join((t1, t2, t3, t4)))
        else:
            dihedral_type = ('-'.join((t4, t3, t2, t1)))
        dihedral_typeids.append(unique_dihedral_types.index(dihedral_type))
        dihedral_groups.append((dihedral.atom1.idx, dihedral.atom2.idx,
                                dihedral.atom3.idx, dihedral.atom4.idx))


    return n_dihedrals, unique_dihedral_types, dihedral_typeids, dihedral_groups
