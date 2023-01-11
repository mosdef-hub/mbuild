"""HOOMD snapshot format."""
import operator
from collections import namedtuple

import numpy as np
import packaging.version
import parmed as pmd

from mbuild.box import Box
from mbuild.compound import Compound, Particle
from mbuild.utils.geometry import coord_shift
from mbuild.utils.io import import_
from mbuild.utils.sorting import natural_sort

hoomd = import_("hoomd")

__all__ = ["to_hoomdsnapshot", "from_snapshot"]


def _get_hoomd_version():
    if "version" in dir(hoomd):
        return packaging.version.parse(hoomd.version.version)
    else:
        return packaging.version.parse(hoomd.__version__)


def from_snapshot(snapshot, scale=1.0):
    """Convert a Snapshot to a Compound.

    Snapshot can be a hoomd.Snapshot or a gsd.hoomd.Snapshot.

    Parameters
    ----------
    snapshot : hoomd.Snapshot or gsd.hoomd.Snapshot
        Snapshot from which to build the mbuild Compound.
    scale : float, optional, default 1.0
        Value by which to scale the length values

    Returns
    -------
    comp : Compound

    Note
    ----
    GSD and HOOMD snapshots center their boxes on the origin (0,0,0), so the
    compound is shifted by half the box lengths
    """
    comp = Compound()
    bond_array = snapshot.bonds.group
    n_atoms = snapshot.particles.N

    if "SnapshotSystemData_float" in dir(hoomd._hoomd) and isinstance(
        snapshot, hoomd._hoomd.SnapshotSystemData_float
    ):
        # hoomd v2
        box = snapshot.box
        comp.box = Box.from_lengths_tilt_factors(
            lengths=np.array([box.Lx, box.Ly, box.Lz]) * scale,
            tilt_factors=np.array([box.xy, box.xz, box.yz]),
        )
    else:
        # gsd / hoomd v3
        box = np.asarray(snapshot.configuration.box)
        comp.box = Box.from_lengths_tilt_factors(
            lengths=box[:3] * scale, tilt_factors=box[3:]
        )

    # GSD and HOOMD snapshots center their boxes on the origin (0,0,0)
    shift = np.array(comp.box.lengths) / 2
    # Add particles
    for i in range(n_atoms):
        name = snapshot.particles.types[snapshot.particles.typeid[i]]
        xyz = snapshot.particles.position[i] * scale + shift
        charge = snapshot.particles.charge[i]

        atom = Particle(name=name, pos=xyz, charge=charge)
        comp.add(atom, label=str(i))

    # Add bonds
    particle_dict = {idx: p for idx, p in enumerate(comp.particles())}
    for i in range(bond_array.shape[0]):
        atom1 = int(bond_array[i][0])
        atom2 = int(bond_array[i][1])
        comp.add_bond([particle_dict[atom1], particle_dict[atom2]])
    return comp


def to_hoomdsnapshot(
    structure,
    ref_distance=1.0,
    ref_mass=1.0,
    ref_energy=1.0,
    rigid_bodies=None,
    shift_coords=True,
    write_special_pairs=True,
    auto_scale=False,
    parmed_kwargs={},
    hoomd_snapshot=None,
):
    """Convert a Compound or parmed.Structure to hoomd.Snapshot.

    Parameters
    ----------
    structure : parmed.Structure
        ParmEd Structure object
        Reference distance for unit conversion ((Angstrom) / (desired units))
    ref_mass : float, optional, default=1.0
        Reference mass for unit conversion ((Dalton) / (desired units))
    ref_energy : float, optional, default=1.0
        Reference energy for unit conversion ((kcal/mol) / (desired units))
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
    hoomd_snapshot : hoomd.Snapshot, optional, default=None
        Initial snapshot to which to add the ParmEd structure object.
        The box information of the initial snapshot will be overwritten.
        (useful for rigid bodies)

    Returns
    -------
    hoomd_snapshot : hoomd.Snapshot
    ReferenceValues : namedtuple
        Values used in scaling

    Notes
    -----
    This method does not create hoomd forcefield objects and
    the snapshot returned does not store the forcefield parameters.
    See mbuild.formats.hoomd_forcefield.create_hoomd_forcefield()

    About units: This method operates on a Parmed.Structure object
        where the units used differ from those used in mBuild and Foyer
        which may have been used when creating the typed Parmed.Structure.

    The default units used when writing out the HOOMD Snapshot are:
    Distance (Angstrom)
    Mass (Dalton)
    Energy (kcal/mol)

    If you wish to convert this unit system to another, you can use the
    reference parameters (ref_distance, ref_mass, ref_energy).
    The values used here should be expected to convert from the Parmed
    Structure units (above) to your desired units.
    The Parmed.Structure values are divided by the reference values.

    If you wish to used a reduced unit system, set auto_scale = True.
    When auto_scale is True, the reference parameters won't be used.

    Examples
    --------
        To convert the energy units from kcal/mol to kj/mol:
            use ref_energy = 0.2390057 (kcal/kj)

        To convert the distance units from Angstrom to nm:
            use ref_distance = 10 (angstroms/nm)

        To use a reduced unit system, where mass, sigma, and epsilon are
        scaled by the largest value of each:
            use auto_scale = True, ref_distance = ref_energy = ref_mass = 1

    """
    if not isinstance(structure, (Compound, pmd.Structure)):
        raise ValueError(
            "You are trying to create a hoomd.Snapshot from a "
            f"{type(structure)} please pass a Compound or pmd.Structure"
        )
    elif isinstance(structure, Compound):
        structure = structure.to_parmed(**parmed_kwargs)

    hoomd_version = _get_hoomd_version()
    if hoomd_version.major == 2 and not hoomd.context.current:
        hoomd.context.initialize("")

    if auto_scale:
        ref_mass = max([atom.mass for atom in structure.atoms])
        pair_coeffs = list(
            set(
                (atom.type, atom.epsilon, atom.sigma)
                for atom in structure.atoms
            )
        )
        ref_energy = max(pair_coeffs, key=operator.itemgetter(1))[1]
        ref_distance = max(pair_coeffs, key=operator.itemgetter(2))[2]

    ReferenceValues = namedtuple("ref_values", ["distance", "mass", "energy"])
    ref_values = ReferenceValues(ref_distance, ref_mass, ref_energy)

    xyz = np.array([[atom.xx, atom.xy, atom.xz] for atom in structure.atoms])
    if shift_coords:
        xyz = coord_shift(xyz, structure.box[:3])

    # Get box information
    if np.allclose(structure.box[3:6], np.array([90, 90, 90])):
        lx, ly, lz = structure.box[:3] / ref_distance
        xy, xz, yz = 0, 0, 0
    else:
        a, b, c = structure.box[0:3] / ref_distance
        alpha, beta, gamma = np.radians(structure.box[3:6])

        lx = a
        xy = b * np.cos(gamma)
        xz = c * np.cos(beta)
        ly = np.sqrt(b**2 - xy**2)
        yz = (b * c * np.cos(alpha) - xy * xz) / ly
        lz = np.sqrt(c**2 - xz**2 - yz**2)

    (
        n_particles,
        scaled_positions,
        unique_types,
        typeids,
        scaled_mass,
        scaled_charges,
        rigid_bodies,
    ) = _parse_particle_information(
        structure, xyz, ref_distance, ref_mass, ref_energy, rigid_bodies
    )
    (
        n_bonds,
        unique_bond_types,
        bond_typeids,
        bond_groups,
    ) = _parse_bond_information(structure)
    (
        n_angles,
        unique_angle_types,
        angle_typeids,
        angle_groups,
    ) = _parse_angle_information(structure)
    (
        n_dihedrals,
        unique_dihedral_types,
        dihedral_typeids,
        dihedral_groups,
    ) = _parse_dihedral_information(structure)
    (
        n_impropers,
        unique_improper_types,
        improper_typeids,
        improper_groups,
    ) = _parse_improper_information(structure)
    pair_types, pair_typeid, pairs, n_pairs = _parse_pair_information(
        structure)

    if hoomd_snapshot is not None:
        n_init = hoomd_snapshot.particles.N

        if n_init > 0:
            n_particles += n_init
            # shift the typeids
            typeids += len(set(hoomd_snapshot.particles.types))
            hoomd_snapshot.particles.types += unique_types
            # shift bond/angle/dihedral indices
            if n_bonds > 0:
                bond_groups = [
                    tuple(row) for row in np.array(bond_groups) + n_init
                ]
                init_bondtypes = len(hoomd_snapshot.bonds.types)
                hoomd_snapshot.bonds.types += unique_bond_types
            if n_angles > 0:
                angle_groups = [
                    tuple(row) for row in np.array(angle_groups) + n_init
                ]
                init_angletypes = len(hoomd_snapshot.angles.types)
                hoomd_snapshot.angles.types += unique_angle_types
            if n_dihedrals > 0:
                dihedral_groups = [
                    tuple(row) for row in np.array(dihedral_groups) + n_init
                ]
                init_dihedraltypes = len(hoomd_snapshot.dihedrals.types)
                hoomd_snapshot.dihedrals.types += unique_dihedral_types
            if n_impropers > 0:
                improper_groups = [
                    tuple(row) for row in np.array(improper_groups) + n_init
                ]
                init_impropertypes = len(hoomd_snapshot.impropers.types)
                hoomd_snapshot.impropers.types += unique_improper_types
            if n_pairs > 0:
                pairs = [tuple(row) for row in np.array(pairs) + n_init]
                init_pairtypes = len(hoomd_snapshot.pairs.types)
                hoomd_snapshot.pairs.types += pair_types
        else:
            raise RuntimeError(
                "Initial snapshot provided, but it contains no particles"
            )

        if hoomd_version.major == 2:
            hoomd_snapshot.box = hoomd.data.boxdim(
                Lx=lx, Ly=ly, Lz=lz, xy=xy, xz=xz, yz=yz
            )

            # save the box for later use when wrapping coordinates
            box = hoomd_snapshot.box
        elif hoomd_version.major == 3:
            hoomd_snapshot.configuration.box = [lx, ly, lz, xy, xz, yz]
        else:
            raise RuntimeError("Unsupported HOOMD version:",
                               str(hoomd_version))

        init_bonds = hoomd_snapshot.bonds.N
        if init_bonds > 0:
            n_bonds += init_bonds
            bond_typeids = list(np.array(bond_typeids) + init_bondtypes)

        init_angles = hoomd_snapshot.angles.N
        if init_angles > 0:
            n_angles += init_angles
            angle_typeids = list(np.array(angle_typeids) + init_angletypes)

        init_dihedrals = hoomd_snapshot.dihedrals.N
        if init_dihedrals > 0:
            n_dihedrals += init_dihedrals
            dihedral_typeids = list(
                np.array(dihedral_typeids) + init_dihedraltypes
            )

        init_impropers = hoomd_snapshot.impropers.N
        if init_impropers > 0:
            n_impropers += init_impropers
            improper_typeids = list(
                np.array(improper_typeids) + init_impropertypes
            )

        init_pairs = hoomd_snapshot.pairs.N
        if init_pairs > 0:
            n_pairs += init_pairs
            pair_typeid = list(np.array(pair_typeid) + init_pairtypes)

    else:
        n_init = 0
        init_bonds = 0
        init_angles = 0
        init_dihedrals = 0
        init_impropers = 0
        init_pairs = 0

        if hoomd_version.major == 2:
            hoomd_snapshot = hoomd.data.make_snapshot(
                N=n_particles,
                box=hoomd.data.boxdim(
                    Lx=lx, Ly=ly, Lz=lz, xy=xy, xz=xz, yz=yz),
                particle_types=unique_types,
                bond_types=unique_bond_types,
                angle_types=unique_angle_types,
                dihedral_types=unique_dihedral_types,
                improper_types=unique_improper_types,
                pair_types=pair_types,
            )
            box = hoomd.data.boxdim(Lx=lx, Ly=ly, Lz=lz, xy=xy, xz=xz, yz=yz)
        elif hoomd_version.major == 3:
            hoomd_snapshot = hoomd.Snapshot()
            hoomd_snapshot.configuration.box = [lx, ly, lz, xy, xz, yz]
            hoomd_snapshot.particles.types = unique_types
            hoomd_snapshot.bonds.types = unique_bond_types
            hoomd_snapshot.angles.types = unique_angle_types
            hoomd_snapshot.dihedrals.types = unique_dihedral_types
            hoomd_snapshot.impropers.types = unique_improper_types
            hoomd_snapshot.pairs.types = pair_types
            box = hoomd.Box(Lx=lx, Ly=ly, Lz=lz, xy=xy, xz=xz, yz=yz)
        else:
            raise RuntimeError("Unsupported HOOMD version:",
                               str(hoomd_version))

    # wrap particles into the box manually for v2
    if hoomd_version.major == 2:
        scaled_positions = np.stack(
            [box.wrap(xyz)[0] for xyz in scaled_positions]
        )

    def set_size(obj, n):
        if hoomd_version.major == 2:
            obj.resize(n)
        elif hoomd_version.major == 3:
            obj.N = n
        else:
            raise RuntimeError("Unsupported HOOMD version:",
                               str(hoomd_version))

    set_size(hoomd_snapshot.particles, n_particles)
    hoomd_snapshot.particles.position[n_init:] = scaled_positions
    hoomd_snapshot.particles.typeid[n_init:] = typeids
    hoomd_snapshot.particles.mass[n_init:] = scaled_mass
    hoomd_snapshot.particles.charge[n_init:] = scaled_charges
    hoomd_snapshot.particles.body[n_init:] = rigid_bodies

    # wrap the particles into the box using the v3 API
    if hoomd_version.major == 3:
        hoomd_snapshot.wrap()

    if n_bonds > 0:
        set_size(hoomd_snapshot.bonds, n_bonds)
        hoomd_snapshot.bonds.typeid[init_bonds:] = bond_typeids
        hoomd_snapshot.bonds.group[init_bonds:] = bond_groups

    if n_angles > 0:
        set_size(hoomd_snapshot.angles, n_angles)
        hoomd_snapshot.angles.typeid[init_angles:] = angle_typeids
        hoomd_snapshot.angles.group[init_angles:] = np.reshape(
            angle_groups, (-1, 3)
        )

    if n_dihedrals > 0:
        set_size(hoomd_snapshot.dihedrals, n_dihedrals)
        hoomd_snapshot.dihedrals.typeid[init_dihedrals:] = dihedral_typeids
        hoomd_snapshot.dihedrals.group[init_dihedrals:] = np.reshape(
            dihedral_groups, (-1, 4)
        )

    if n_impropers > 0:
        set_size(hoomd_snapshot.impropers, n_impropers)
        hoomd_snapshot.impropers.typeid[init_impropers:] = improper_typeids
        hoomd_snapshot.impropers.group[init_impropers:] = np.reshape(
            improper_groups, (-1, 4)
        )

    if n_pairs > 0:
        set_size(hoomd_snapshot.pairs, n_pairs)
        hoomd_snapshot.pairs.typeid[init_pairs:] = pair_typeid
        hoomd_snapshot.pairs.group[init_pairs:] = np.reshape(pairs, (-1, 2))

    return hoomd_snapshot, ref_values


def _parse_particle_information(
    structure, xyz, ref_distance, ref_mass, ref_energy, rigid_bodies
):
    """Write out the particle information."""
    n_particles = len(structure.atoms)
    scaled_positions = xyz / ref_distance

    types = [
        atom.name if atom.type == "" else atom.type for atom in structure.atoms
    ]

    unique_types = list(set(types))
    unique_types.sort(key=natural_sort)

    typeids = np.array([unique_types.index(t) for t in types])

    masses = np.array([atom.mass for atom in structure.atoms])
    masses[masses == 0] = 1.0
    scaled_mass = masses / ref_mass

    charges = np.array([atom.charge for atom in structure.atoms])
    e0 = 2.396452e-04
    """
    Permittivity of free space = 2.396452e-04 e^2/((kcal/mol) Angstrom),
    where e is the elementary charge
    """
    charge_factor = (4.0 * np.pi * e0 * ref_distance * ref_energy) ** 0.5
    scaled_charges = charges / charge_factor

    if rigid_bodies:
        rigid_bodies = [-1 if body is None else body for body in rigid_bodies]
    else:
        rigid_bodies = [-1 for _ in structure.atoms]

    return (
        n_particles,
        scaled_positions,
        unique_types,
        typeids,
        scaled_mass,
        scaled_charges,
        rigid_bodies,
    )


def _parse_pair_information(structure):
    pair_types = []
    pair_typeid = []
    pairs = []
    for ai in structure.atoms:
        for aj in ai.dihedral_partners:
            # make sure we don't double add
            if ai.idx > aj.idx:
                ps = "-".join(sorted([ai.type, aj.type]))
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
        if t1 == "" or t2 == "":
            t1, t2 = bond.atom1.name, bond.atom2.name
        t1, t2 = sorted([t1, t2], key=natural_sort)
        try:
            bond_type = "-".join((t1, t2))
        except AttributeError:  # no forcefield applied, bond.type is None
            bond_type = ("-".join((t1, t2)), 0.0, 0.0)
        unique_bond_types.add(bond_type)
    unique_bond_types = sorted(list(unique_bond_types), key=natural_sort)

    bond_typeids = []
    bond_groups = []
    for bond in structure.bonds:
        t1, t2 = bond.atom1.type, bond.atom2.type
        if t1 == "" or t2 == "":
            t1, t2 = bond.atom1.name, bond.atom2.name
        t1, t2 = sorted([t1, t2], key=natural_sort)
        try:
            bond_type = "-".join((t1, t2))
        except AttributeError:  # no forcefield applied, bond.type is None
            bond_type = ("-".join((t1, t2)), 0.0, 0.0)
        bond_typeids.append(unique_bond_types.index(bond_type))
        bond_groups.append((bond.atom1.idx, bond.atom2.idx))

    return n_bonds, unique_bond_types, bond_typeids, bond_groups


def _parse_angle_information(structure):
    n_angles = len(structure.angles)

    unique_angle_types = set()
    for angle in structure.angles:
        t1, t2, t3 = angle.atom1.type, angle.atom2.type, angle.atom3.type
        t1, t3 = sorted([t1, t3], key=natural_sort)
        angle_type = "-".join((t1, t2, t3))
        unique_angle_types.add(angle_type)
    unique_angle_types = sorted(list(unique_angle_types), key=natural_sort)

    angle_typeids = []
    angle_groups = []
    for angle in structure.angles:
        t1, t2, t3 = angle.atom1.type, angle.atom2.type, angle.atom3.type
        t1, t3 = sorted([t1, t3], key=natural_sort)
        angle_type = "-".join((t1, t2, t3))
        angle_typeids.append(unique_angle_types.index(angle_type))
        angle_groups.append(
            (angle.atom1.idx, angle.atom2.idx, angle.atom3.idx))

    return n_angles, unique_angle_types, angle_typeids, angle_groups


def _parse_dihedral_information(structure):
    n_dihedrals = len(structure.rb_torsions + structure.dihedrals)

    unique_dihedral_types = set()
    for dihedral in structure.rb_torsions + structure.dihedrals:
        t1, t2 = dihedral.atom1.type, dihedral.atom2.type
        t3, t4 = dihedral.atom3.type, dihedral.atom4.type
        if [t2, t3] == sorted([t2, t3], key=natural_sort):
            dihedral_type = "-".join((t1, t2, t3, t4))
        else:
            dihedral_type = "-".join((t4, t3, t2, t1))
        unique_dihedral_types.add(dihedral_type)
    unique_dihedral_types = sorted(
        list(unique_dihedral_types), key=natural_sort
    )

    dihedral_typeids = []
    dihedral_groups = []
    for dihedral in structure.rb_torsions + structure.dihedrals:
        t1, t2 = dihedral.atom1.type, dihedral.atom2.type
        t3, t4 = dihedral.atom3.type, dihedral.atom4.type
        if [t2, t3] == sorted([t2, t3], key=natural_sort):
            dihedral_type = "-".join((t1, t2, t3, t4))
        else:
            dihedral_type = "-".join((t4, t3, t2, t1))
        dihedral_typeids.append(unique_dihedral_types.index(dihedral_type))
        dihedral_groups.append(
            (
                dihedral.atom1.idx,
                dihedral.atom2.idx,
                dihedral.atom3.idx,
                dihedral.atom4.idx,
            )
        )

    return n_dihedrals, unique_dihedral_types, dihedral_typeids, dihedral_groups


def _parse_improper_information(structure):
    n_dihedrals = len(structure.impropers)

    unique_dihedral_types = set()
    for dihedral in structure.impropers:
        t1, t2 = dihedral.atom1.type, dihedral.atom2.type
        t3, t4 = dihedral.atom3.type, dihedral.atom4.type
        if [t2, t3] == sorted([t2, t3], key=natural_sort):
            dihedral_type = "-".join((t1, t2, t3, t4))
        else:
            dihedral_type = "-".join((t4, t3, t2, t1))
        unique_dihedral_types.add(dihedral_type)
    unique_dihedral_types = sorted(
        list(unique_dihedral_types), key=natural_sort
    )

    dihedral_typeids = []
    dihedral_groups = []
    for dihedral in structure.impropers:
        t1, t2 = dihedral.atom1.type, dihedral.atom2.type
        t3, t4 = dihedral.atom3.type, dihedral.atom4.type
        if [t2, t3] == sorted([t2, t3], key=natural_sort):
            dihedral_type = "-".join((t1, t2, t3, t4))
        else:
            dihedral_type = "-".join((t4, t3, t2, t1))
        dihedral_typeids.append(unique_dihedral_types.index(dihedral_type))
        dihedral_groups.append(
            (
                dihedral.atom1.idx,
                dihedral.atom2.idx,
                dihedral.atom3.idx,
                dihedral.atom4.idx,
            )
        )

    return n_dihedrals, unique_dihedral_types, dihedral_typeids, dihedral_groups
